
# from .models.utils import modified_exponential
import math as m
import vice



# --------------- RADIAL GAS FLOWS --------------- #
RADIAL_GAS_FLOWS = None # None turns them off
RADIAL_GAS_FLOW_ONSET = 1 # Gyr -- radial flow starts 1 Gyr in

# used when RADIAL_GAS_FLOWS = "constant"
RADIAL_GAS_FLOW_SPEED = -1.5 # km/s
# def RADIAL_GAS_FLOW_SPEED(time):
# 	return -10 * np.exp(-time / 3)

# used when RADIAL_GAS_FLOWS = "angular_momentum_dilution"
# RADIAL_GAS_FLOW_BETA_PHI_IN = 0.7
RADIAL_GAS_FLOW_BETA_PHI_IN = 0.3
# def RADIAL_GAS_FLOW_BETA_PHI_IN(r, t):
	# return 0.3 + 0.4 * (1 - m.exp(-t / 2))
RADIAL_GAS_FLOW_BETA_PHI_OUT = 0


# --------------- ACCRETION METALLICITY TIME-DEP --------------- #
CGM_FINAL_METALLICITY = -float("inf") # -inf for zero metallicity accretion
CGM_METALLICITY_GROWTH_TIMESCALE = 3


# --------------- YIELDS --------------- #
YIELDSOLAR = 1
FE_CC_FRAC = 0.35

vice.yields.ccsne.settings["o"] = YIELDSOLAR * vice.solar_z["o"]
vice.yields.sneia.settings["o"] = 0
vice.yields.ccsne.settings["fe"] = FE_CC_FRAC * YIELDSOLAR * vice.solar_z["fe"]
vice.yields.sneia.settings["fe"] = (1 - FE_CC_FRAC) * YIELDSOLAR * vice.solar_z["fe"]

