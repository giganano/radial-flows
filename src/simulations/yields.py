
import vice

YIELDSOLAR = 1
FE_CC_FRAC = 0.35

vice.yields.ccsne.settings["o"] = YIELDSOLAR * vice.solar_z["o"]
vice.yields.sneia.settings["o"] = 0
vice.yields.ccsne.settings["fe"] = FE_CC_FRAC * YIELDSOLAR * vice.solar_z["fe"]
vice.yields.sneia.settings["fe"] = (1 - FE_CC_FRAC) * YIELDSOLAR * vice.solar_z["fe"]

