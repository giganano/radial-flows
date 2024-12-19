
import vice



# --------------- RADIAL GAS FLOWS --------------- #
RADIAL_GAS_FLOWS = "river" # Turn them off
RADIAL_GAS_FLOW_ONSET = 1 # Gyr -- radial flow starts 1 Gyr in

# used when RADIAL_GAS_FLOWS = "constant"
RADIAL_GAS_FLOW_SPEED = -1 # km/s

# used when RADIAL_GAS_FLOWS = "angular_momentum_dilution"
RADIAL_GAS_FLOW_BETA_PHI_IN = 0.7
RADIAL_GAS_FLOW_BETA_PHI_OUT = 0




# --------------- YIELDS --------------- #
YIELDSOLAR = 1
FE_CC_FRAC = 0.35

vice.yields.ccsne.settings["o"] = YIELDSOLAR * vice.solar_z["o"]
vice.yields.sneia.settings["o"] = 0
vice.yields.ccsne.settings["fe"] = FE_CC_FRAC * YIELDSOLAR * vice.solar_z["fe"]
vice.yields.sneia.settings["fe"] = (1 - FE_CC_FRAC) * YIELDSOLAR * vice.solar_z["fe"]

