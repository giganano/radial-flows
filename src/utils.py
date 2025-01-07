
import numpy as np

def get_velocity_profile(output, lookback):
	raw = np.genfromtxt("%s_gasvelocities.out" % (output.name))
	time = output.zones["zone0"].history["time"][-1] - lookback
	diff = [abs(_ - time) for _ in output.zones["zone0"].history["time"]]
	idx = diff.index(min(diff))
	time = output.zones["zone0"].history["time"][idx]
	radii = []
	vgas = []
	for i in range(len(raw)):
		if raw[i][0] == time:
			radii.append(raw[i][1])
			vgas.append(raw[i][2])
		else: pass
	return [radii, vgas]

def mu(output, lookback, zone_width = 0.1):
	radii, vgas = get_velocity_profile(output, lookback)
	diff = [abs(_ - lookback) for _ in output.zones["zone0"].history["lookback"]]
	idx = diff.index(min(diff))
	mu_gas = []
	mu_oxygen = []
	for i in range(len(radii) - 1):
		zone = output.zones["zone%d" % (i)]
		neighbor = output.zones["zone%d" % (i + 1)]
		if radii[i + 1] >= 15.5:
			mu_gas.append(float("nan"))
			mu_oxygen.append(float("nan"))
		else:
			tau_star = zone.history["mgas"][idx] / zone.history["sfr"][idx] * 1.e-9
			mu = (neighbor.history["mgas"][idx] - zone.history["mgas"][idx]) / (
				zone.history["mgas"][idx] * zone_width)
			mu += (vgas[i + 1] - vgas[i]) / (vgas[i] * zone_width)
			mu *= -tau_star * vgas[i]
			mu_gas.append(mu)
			mu -= tau_star * vgas[i] * (
				neighbor.history["z(o)"][idx] - zone.history["z(o)"][idx]) / (
				zone.history["z(o)"][idx] * zone_width)
			mu_oxygen.append(mu)
	return [radii, mu_gas, mu_oxygen]
