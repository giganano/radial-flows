
from scipy import interpolate
from .inputs import YIELDSOLAR
import numpy as np
import vice
import os

# adjust primordial abundances
D_MASS_RATIO = 2.014 / 1.007
vice.primordial["he"] = 0.24719 # helium-4
vice.primordial["au"] = 2.341e-5 # helium-3
# deuterium (Ryan Cooke et al. 2014)
vice.primordial["ag"] = (1 - vice.primordial["he"] - vice.primordial["au"]) * D_MASS_RATIO * 2.53e-5
vice.elements.nonmetals("he", "au", "ag")
vice.elements.destroyed_by_stars("ag")

yieldsdir = "%s/" % (os.path.abspath(os.path.dirname(__file__)))

class amplified_ccsn(interpolate.interp1d):

    def __init__(self, *args, prefactor = 1, **kwargs):
        super().__init__(*args, **kwargs)
        self.prefactor = prefactor

    def __call__(self, z):
        result = self.prefactor * super().__call__([z])[0]
        return result

# Setup AGB yields
mass3, Zmet3, yldpm3 = np.loadtxt(yieldsdir + "Lagarde2011_3He_yields.csv",
    delimiter=",", unpack=True, usecols=(0,1,3))
mass4, Zmet4, yldpm4 = np.loadtxt(yieldsdir + "Lagarde2011_4He_yields.csv",
    delimiter=",", unpack=True, usecols=(0,1,3))
sh = (4,9,)
metval = Zmet3.reshape(sh)[:,0]
# Add a value for 8 Msun
massval = np.append(mass3.reshape(sh)[0,:], 8.0)
yld3 = np.append(yldpm3.reshape(sh).T, -1.75E-5*np.ones((1,metval.size)), axis=0)
yld4 = np.append(yldpm4.reshape(sh).T, yldpm4.reshape(sh).T[-1,:].reshape((1,metval.size)), axis=0)
yld3spl_agb = interpolate.RectBivariateSpline(massval, metval, yld3, kx=1, ky=1)
yld4spl_agb = interpolate.RectBivariateSpline(massval, metval, yld4, kx=1, ky=1)
# Don't add a value for 8 Msun
# massval = mass3.reshape(sh)[0,:]
# yld3spl_agb = interpolate.RectBivariateSpline(massval, metval, yldpm3.reshape(sh).T, kx=1, ky=1)
# yld4spl_agb = interpolate.RectBivariateSpline(massval, metval, yldpm4.reshape(sh).T, kx=1, ky=1)

# Setup CCSNe yields
data3 = np.load(yieldsdir + "LC2018_3He_yields_Kroupa.npy")
data4 = np.load(yieldsdir + "LC2018_4He_yields_Kroupa.npy")
# data3 = np.load(yieldsdir + "CL2004_3He_yields_Kroupa.npy")
# data4 = np.load(yieldsdir + "CL2004_4He_yields_Kroupa.npy")
#data3 = np.load(yieldsdir + "CL2004_3He_yields_Scalo86.npy")
#data4 = np.load(yieldsdir + "CL2004_4He_yields_Scalo86.npy")
yld3spl_ccsne = interpolate.interp1d(data3[:,0], data3[:,1], kind='linear', bounds_error=False, fill_value='extrapolate')
yld4spl_ccsne = amplified_ccsn(data4[:,0], data4[:,1], kind='linear', bounds_error=False, fill_value='extrapolate')


def custom_agb_yield_3He(mass, z):
    # Mass and metallicity dependent model
    return yld3spl_agb(mass, z)[0][0]


def custom_agb_yield_4He(mass, z):
    # Mass and metallicity dependent model
    return yld4spl_agb(mass, z)[0][0]


def custom_ccsne_yield_3He(zval):
    # Metallicity dependent model
    return yld3spl_ccsne([zval])[0]


# def custom_ccsne_yield_4He(zval):
#     # Metallicity dependent model
#     return YIELDSOLAR * yld4spl_ccsne([zval])[0]


def zero_agb_yield(mass, z):
    return 0


# helium-4
vice.yields.ccsne.settings['he'] = yld4spl_ccsne
vice.yields.agb.settings['he'] = custom_agb_yield_4He


# helium-3
vice.yields.ccsne.settings['au'] = custom_ccsne_yield_3He
vice.yields.sneia.settings['au'] = 0
vice.yields.agb.settings['au'] = custom_agb_yield_3He

# deuterium
vice.yields.ccsne.settings['ag'] = 0
vice.yields.sneia.settings['ag'] = 0
vice.yields.agb.settings['ag'] = zero_agb_yield


# Set a metallicity dependent mass-lifetime relation (MLR)
#vice.mlr.setting = "hpt2000"

# adjust normalization of helium-4 yield to Weller et al. (2025)
current_he4_ccsn = vice.yields.ccsne.settings['he']
current_he4_agb = vice.yields.agb.settings['he']
vice.yields.ccsne.settings['he'] = 0
agbmass, times = vice.single_stellar_population('he', mstar = 1.0e6, time = 13.2)
vice.yields.ccsne.settings['he'] = current_he4_ccsn
vice.yields.agb.settings['he'] = zero_agb_yield
ccsnmass, _ = vice.single_stellar_population('he', mstar = 1.0e6, time = 13.2)
vice.yields.ccsne.settings['he'] = current_he4_ccsn
vice.yields.agb.settings['he'] = current_he4_agb
target = 0.022 * 0.0073 / vice.solar_z['o']
target *= vice.yields.ccsne.settings['o'] / 0.0071
target *= 1.0e6
vice.yields.ccsne.settings['he'].prefactor *= (target - agbmass[-1]) / ccsnmass[-1]
# print(vice.yields.ccsne.settings['he'].prefactor)







