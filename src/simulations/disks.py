r"""
The diskmodel objects employed in the Johnson et al. (2021) study.
"""

try:
	ModuleNotFoundError
except NameError:
	ModuleNotFoundError = ImportError
try:
	import vice
except (ModuleNotFoundError, ImportError):
	raise ModuleNotFoundError("Could not import VICE.")
if vice.version[:2] < (1, 2):
	raise RuntimeError("""VICE version >= 1.2.0 is required to produce \
Johnson et al. (2021) figures. Current: %s""" % (vice.__version__))
else: pass
from vice.yields.presets import JW20
from vice.toolkit import hydrodisk
vice.yields.sneia.settings['fe'] *= 10**0.1
from .._globals import END_TIME, MAX_SF_RADIUS, ZONE_WIDTH
from . import gasflows
from . import outflows
from . import migration
from . import models
from . import inputs
from . import sfe
from .models.utils import get_bin_number, interpolate, modified_exponential
from .models.gradient import gradient
# import warnings
import math as m
import sys



class diskmodel(vice.milkyway):

	r"""
	A milkyway object tuned to the Johnson et al. (2021) models specifically.

	Parameters
	----------
	zone_width : ``float`` [default : 0.1]
		The width of each annulus in kpc.
	name : ``str`` [default : "diskmodel"]
		The name of the model; the output will be stored in a directory under
		this name with a ".vice" extension.
	spec : ``str`` [default : "static"]
		A keyword denoting the time-dependence of the star formation history.
		Allowed values:

		- "static"
		- "insideout"
		- "lateburst"
		- "outerburst"

	verbose : ``bool`` [default : True]
		Whether or not the run the models with verbose output.
	migration_mode : ``str`` [default : "diffusion"]
		A keyword denoting the time-dependence of stellar migration.
		Allowed values:

		- "diffusion"
		- "linear"
		- "sudden"
		- "post-process"

	kwargs : varying types
		Other keyword arguments to pass ``vice.milkyway``.

	Attributes and functionality are inherited from ``vice.milkyway``.
	"""

	def __init__(self, zone_width = 0.1, timestep_size = 0.01,
		name = "diskmodel", spec = "static", verbose = True,
		migration_mode = "diffusion", **kwargs):
		super().__init__(zone_width = zone_width, name = name,
			verbose = verbose, **kwargs)
		if self.zone_width <= 0.2 and self.dt <= 0.02 and self.n_stars >= 6:
			Nstars = 3102519
		else:
			Nstars = 2 * int(MAX_SF_RADIUS / zone_width * END_TIME / self.dt *
				self.n_stars)
		# self.migration.stars = migration.diskmigration(self.annuli,
		# 	N = Nstars, mode = migration_mode,
		# 	filename = "%s_analogdata.out" % (name))
		self.migration.stars = migration.gaussian_migration(self.annuli,
			zone_width = zone_width,
			filename = "%s_analogdata.out" % (self.name),
			post_process = self.simple)
		# self.migration.stars = migration.nomigration(self.annuli,
		# 	zone_width = zone_width,
		# 	filename = "%s_analogdata.out" % (self.name),
		# 	post_process = self.simple)
		self.evolution = star_formation_history(spec = spec,
			zone_width = zone_width)
		self.mode = "sfr"
		self.dt = timestep_size

		for i in range(self.n_zones):
			if inputs.OUTFLOWS in ["empirical_calib", "J25", "rc25_constant"]:
				cls = {
					"empirical_calib": outflows.empirical_calib,
					"J25": outflows.J25,
					"rc25_constant": outflows.rc25_constant,
				}[inputs.OUTFLOWS]
				kwargs = {
					# "zone_width": zone_width,
					"timestep": self.zones[i].dt
				}
				if i: kwargs["evol"] = self.zones[0].eta.evol
				self.zones[i].eta = cls(self,
					i * zone_width + 1.e-6, **kwargs)
			elif inputs.OUTFLOWS == "constant_t_and_r":
				self.zones[i].eta = outflows.constant_t_and_r(
					inputs.OUTFLOWS_CONST_ETA)
			elif inputs.OUTFLOWS is None:
				self.zones[i].eta = 0
			else:
				raise ValueError("Bad outflow setting in input file.")

		for i in range(self.n_zones):
			area = m.pi * ZONE_WIDTH**2 * ((i + 1)**2 - i**2)
			self.zones[i].tau_star = sfe.sfe(area, mode = "sfr")
			# self.zones[i].tau_star = lambda r, t: 3

		# setup radial gas flow
		if inputs.RADIAL_GAS_FLOWS is not None:
			kwargs = {
				"onset": inputs.RADIAL_GAS_FLOW_ONSET,
				"dr": zone_width,
				"dt": self.dt,
				"outfilename": "%s_gasvelocities.out" % (self.name)
			}
			callkwargs = {}
			if inputs.RADIAL_GAS_FLOWS == "constant":
				engine = gasflows.constant(inputs.RADIAL_GAS_FLOW_SPEED,
					**kwargs)
			elif inputs.RADIAL_GAS_FLOWS == "oscillatory":
				engine = gasflows.oscillatory(
					inputs.RADIAL_GAS_FLOW_MEAN,
					inputs.RADIAL_GAS_FLOW_AMPLITUDE,
					inputs.RADIAL_GAS_FLOW_PERIOD,
					**kwargs)
			elif inputs.RADIAL_GAS_FLOWS == "linear":
				engine = gasflows.linear(dvdr = inputs.RADIAL_GAS_FLOW_DVDR,
					**kwargs)
			elif inputs.RADIAL_GAS_FLOWS == "angular_momentum_dilution":
				engine = gasflows.angular_momentum_dilution(self,
					beta_phi_in = inputs.RADIAL_GAS_FLOW_BETA_PHI_IN,
					beta_phi_out = inputs.RADIAL_GAS_FLOW_BETA_PHI_OUT,
					**kwargs)
				# callkwargs["recycling"] = 0.4
			elif inputs.RADIAL_GAS_FLOWS == "potential_well_deepening":
				engine = gasflows.potential_well_deepening(self,
					gamma = inputs.RADIAL_GAS_FLOW_PWDGAMMA, **kwargs)
			elif inputs.RADIAL_GAS_FLOWS == "river":
				engine = gasflows.river(self, **kwargs)
				callkwargs["recycling"] = 0.4
			else:
				raise ValueError(
					"Unrecognized radial gas flow setting: %s" % (
						inputs.RADIAL_GAS_FLOWS))
			engine.setup(self, **callkwargs)
		else:
			pass

		if not m.isinf(inputs.CGM_FINAL_METALLICITY):
			for i in range(self.n_zones):
				self.zones[i].Zin = {}
				for elem in self.zones[i].elements:
					kwargs = {
						"norm": vice.solar_z[elem] * 10**inputs.CGM_FINAL_METALLICITY,
						"rise": inputs.CGM_METALLICITY_GROWTH_TIMESCALE,
						"timescale": float("inf")
					}
					self.zones[i].Zin[elem] = modified_exponential(**kwargs)
		else: pass


		# setup radial gas flow
# 		vgas_engine = gasflows.river(self,
# 			outfilename = "%s_gasvelocities.out" % (self.name))
# 		vgas_alltimes = []
# 		times = [self.dt * i for i in range(int(END_TIME / self.dt) + 10)]
# 		for i in range(len(times)):
# 			if i > MIGRATION_PAUSE / self.dt:
# 				radii, vgas = vgas_engine(i * self.dt)
# 				vgas_alltimes.append(vgas)
# 			else:
# 				radii = [zone_width * i for i in range(self.n_zones)]
# 				vgas = len(radii) * [0.]
# 				vgas_alltimes.append(vgas)
# 		matrix_elements_inward = []
# 		matrix_elements_outward = []
# 		for i in range(self.n_zones):
# 			yvals_inward = []
# 			yvals_outward = []
# 			vgas = [row[i] for row in vgas_alltimes]
# 			for j in range(len(times)):
# 				if j > MIGRATION_PAUSE / self.dt:
# 					radius = i * zone_width
# 					if vgas[j] > 0: # outward flow
# 						numerator = 2 * (radius + zone_width) * vgas[j] * 0.01
# 						numerator -= vgas[j]**2 * 0.01**2
# 					else: # inward flow
# 						numerator = vgas[j]**2 * 0.01**2
# 						numerator -= 2 * radius * vgas[j] * 0.01
# 						# numerator -= 2 * radius * zone_width + zone_width**2
# 					denominator = 2 * radius * zone_width + zone_width**2
# 					areafrac = numerator / denominator
# 					if areafrac * self.dt / 0.01 > 1:
# 						warnings.warn("""\
# Area fraction larger than 1. Consider comparing results with different \
# timestep sizes to assess the impact of numerical artifacts.""")
# 						areafrac = 0.01 / self.dt - 1.e-9
# 					elif areafrac < 0:
# 						areafrac = 1.e-9
# 					else: pass
# 					if vgas[j] > 0:
# 						yvals_outward.append(areafrac)
# 						yvals_inward.append(1.e-10)
# 					else:
# 						yvals_outward.append(1.e-10)
# 						yvals_inward.append(areafrac)
# 				else:
# 					yvals_outward.append(1.e-10)
# 					yvals_inward.append(1.e-10)
# 			matrix_elements_outward.append(
# 				gasflows.driver(times, yvals_outward, dt = self.dt))
# 			matrix_elements_inward.append(
# 				gasflows.driver(times, yvals_inward, dt = self.dt))
# 		for i in range(self.n_zones):
# 			for j in range(self.n_zones):
# 				if i - 1 == j: # inward flows
# 					self.migration.gas[i][j] = matrix_elements_inward[i]
# 				elif i + 1 == j: # outward flows
# 					self.migration.gas[i][j] = matrix_elements_outward[i]
# 				else:
# 					self.migration.gas[i][j] = 0



	def run(self, *args, **kwargs):
		out = super().run(*args, **kwargs)
		self.migration.stars.close_file()
		return out

	@classmethod
	def from_config(cls, config, **kwargs):
		r"""
		Obtain a ``diskmodel`` object with the parameters encoded into a
		``config`` object.

		**Signature**: diskmodel.from_config(config, **kwargs)

		Parameters
		----------
		config : ``config``
			The ``config`` object with the parameters encoded as attributes.
			See src/simulations/config.py.
		**kwargs : varying types
			Additional keyword arguments to pass to ``diskmodel.__init__``.

		Returns
		-------
		model : ``diskmodel``
			The ``diskmodel`` object with the proper settings.
		"""
		model = cls(zone_width = config.zone_width,
			timestep_size = config.timestep_size, **kwargs)
		model.n_stars = config.star_particle_density
		model.bins = config.bins
		model.elements = config.elements
		return model


class star_formation_history:

	r"""
	The star formation history (SFH) of the model galaxy. This object will be
	used as the ``evolution`` attribute of the ``diskmodel``.

	Parameters
	----------
	spec : ``str`` [default : "static"]
		A keyword denoting the time-dependence of the SFH.
	zone_width : ``float`` [default : 0.1]
		The width of each annulus in kpc.

	Calling
	-------
	- Parameters

		radius : ``float``
			Galactocentric radius in kpc.
		time : ``float``
			Simulation time in Gyr.
	"""

	def __init__(self, spec = "static", zone_width = 0.1):
		self._radii = []
		self._evol = []
		i = 0
		max_radius = 20 # kpc, defined by ``vice.milkyway`` object.
		while (i + 1) * zone_width < max_radius:
			self._radii.append((i + 0.5) * zone_width)
			self._evol.append({
					"oscil":		models.insideout_oscil,
					"static": 		models.static,
					"insideout": 	models.insideout,
					"lateburst": 	models.lateburst,
					"outerburst": 	models.outerburst
				}[spec.lower()]((i + 0.5) * zone_width))
			i += 1

	def __call__(self, radius, time):
		# The milkyway object will always call this with a radius in the
		# self._radii array, but this ensures a continuous function of radius
		if radius > MAX_SF_RADIUS:
			return 0
		else:
			idx = get_bin_number(self._radii, radius)
			if idx != -1:
				return gradient(radius) * interpolate(self._radii[idx],
					self._evol[idx](time), self._radii[idx + 1],
					self._evol[idx + 1](time), radius)
			else:
				return gradient(radius) * interpolate(self._radii[-2],
					self._evol[-2](time), self._radii[-1], self._evol[-1](time),
					radius)

