r"""
Handles radial gas flows in these models.
"""

from .._globals import MAX_SF_RADIUS, END_TIME
from .models.utils import get_bin_number, sinusoid, logistic
from .outflows import evoldata
from . import inputs
from vice.toolkit.interpolation import interp_scheme_1d
from vice.milkyway.milkyway import _MAX_RADIUS_ as MAX_RADIUS # 20 kpc
from scipy.integrate import solve_ivp
import numpy as np
import warnings
import vice
from vice import ScienceWarning


class driver(interp_scheme_1d):

	def __init__(self, *args, dt = 0.01, **kwargs):
		super().__init__(*args, **kwargs)
		self.dt = dt
	
	def __call__(self, x):
		# test = super().__call__(x) * self.dt / 0.01
		# if test < 0:
		# 	return 0
		# elif test > 1:
		# 	return 1
		# else:
		# 	return test

		test = super().__call__(x)
		if test < 0:
			return 0
		else:
			return min(test, 0.01 / self.dt - 1.e-9)
		# elif test * self.dt / 0.01 > 1 - 1.e-9:
		# 	return 0.01 / self.dt - 1.e-9
		# else:
		# 	return test



class base:

	NORMTIME = 0.01 # Gyr

	def __init__(self, onset = 1, dr = 0.1, dt = 0.01,
		outfilename = "gasvelocities.out"):
		if outfilename is not None:
			self.outfile = open(outfilename, "w")
			self.outfile.write("# Time [Gyr]    ")
			self.outfile.write("Radius [kpc]    ")
			self.outfile.write("ISM radial velocity [kpc/Gyr]\n")
		else:
			self.outfile = None
		self.onset = onset
		self.dr = dr
		self.dt = dt
		# if bar:
		# 	self.bar = barflows(inputs.BAR_SPEED,
		# 		inneredge = inputs.BAR_INNER,
		# 		outeredge = inputs.BAR_OUTER,
		# 		innerscale = inputs.BAR_INNER_SCALE,
		# 		outerscale = inputs.BAR_OUTER_SCALE)
		# else:
		# 	self.bar = None


	def __enter__(self):
		return self


	def __exit__(self, exc_type, exc_value, exc_tb):
		self.outfile.close()
		return exc_value is None


	def write(self, time, radii, velocities):
		if self.outfile is not None:
			for i in range(len(radii)):
				self.outfile.write("%.5e\t%.5e\t%.5e\n" % (
					time, radii[i], velocities[i]))
		else: pass


	def setup(self, mw_model, **kwargs):
		vgas_alltimes = []
		n_zones = int(MAX_RADIUS / self.dr)
		# n_zones = len(mw_model.annuli)
		times = [self.dt * i for i in range(int(END_TIME / self.dt) + 10)]
		for i in range(len(times)):
			if i > self.onset / self.dt:
				radii, vgas = self.__call__(i * self.dt, **kwargs)
				vgas_alltimes.append(vgas)
			else:
				radii = [self.dr * i for i in range(n_zones)]
				vgas = len(radii) * [0.]
				vgas_alltimes.append(vgas)
		matrix_elements_inward = []
		matrix_elements_outward = []
		for i in range(n_zones):
			areafracs_inward = []
			areafracs_outward = []
			vgas = [row[i] for row in vgas_alltimes]
			for j in range(len(times)):
				if j >= self.onset / self.dt:
					radius = i * self.dr
					if vgas[j] > 0: # outward flow
						numerator = 2 * (radius + 
							self.dr) * vgas[j] * self.NORMTIME
						numerator -= vgas[j]**2 * self.NORMTIME**2
					else: # inward flow
						numerator = vgas[j]**2 * self.NORMTIME**2
						numerator -= 2 * radius * vgas[j] * self.NORMTIME
					denominator = 2 * radius * self.dr + self.dr**2
					areafrac = numerator / denominator
					if areafrac * self.dt / self.NORMTIME > 1 - 1.e-9:
						warnings.warn("""\
Area fraction larger than 1. Consider comparing results with different \
timestep sizes to assess the impact of numerical artifacts.""", ScienceWarning)
						areafrac = 1 - 1.e-9
					elif areafrac * self.dt / self.NORMTIME < 1.e-9:
						areafrac = 1.e-9
					else: pass
					if vgas[j] > 0:
						areafracs_outward.append(areafrac)
						areafracs_inward.append(1.e-10)
					else:
						areafracs_outward.append(1.e-10)
						areafracs_inward.append(areafrac)
				else:
					areafracs_outward.append(1.e-10)
					areafracs_inward.append(1.e-10)
			matrix_elements_outward.append(
				driver(times, areafracs_outward, dt = self.dt))
			matrix_elements_inward.append(
				driver(times, areafracs_inward, dt = self.dt))
		for i in range(n_zones):
			for j in range(n_zones):
				if i - 1 == j: # inward flows
					mw_model.migration.gas[i][j] = matrix_elements_inward[i]
				elif i + 1 == j: # outward flows
					mw_model.migration.gas[i][j] = matrix_elements_outward[i]
				else:
					mw_model.migration.gas[i][j] = 0




class constant(base):

	def __init__(self, speed, onset = 1, dr = 0.1, dt = 0.01,
		outfilename = "gasvelocities.out"):
		super().__init__(onset = onset, dr = dr, dt = dt,
			outfilename = outfilename)
		self.speed = speed


	def __call__(self, time, **kwargs):
		radii = [self.dr * i for i in range(int(MAX_RADIUS / self.dr))]
		if callable(self.speed):
			# it's a constant in radius, but not necessarily in time
			speed = self.speed(time, **kwargs)
		else:
			speed = self.speed
		vgas = len(radii) * [speed]
		self.write(time, radii, vgas)
		return [radii, vgas]


class steady(base, logistic):

	def __init__(self, maximum, midpoint = 3, scale = 0.5,
		onset = 1, dr = 0.1, dt = 0.01, outfilename = "gasvelocities.out"):
		base.__init__(self, onset = onset, dr = dr, dt = dt,
			outfilename = outfilename)
		logistic.__init__(self, maximum = maximum, midpoint = midpoint,
			scale = scale)

	def __call__(self, time, **kwargs):
		radii = [self.dr * i for i in range(int(MAX_RADIUS / self.dr))]
		vgas = [logistic.__call__(self, r) for r in radii]
		self.write(time, radii, vgas)
		return [radii, vgas]


class limexp(constant):

	def __call__(self, time, **kwargs):
		radii = [self.dr * i for i in range(int(MAX_RADIUS / self.dr))]
		if callable(self.speed):
			speed = self.speed(time, **kwargs)
		else:
			speed = self.speed
		vgas = []
		for i in range(len(radii)):
			vgas.append(speed * (1 - np.exp(-radii[i] / 3)))
		self.write(time, radii, vgas)
		return [radii, vgas]


# class steady(base):

# 	def __init__(self, maximum, scale = 6,
# 		onset = 1, dr = 0.1, dt = 0.01, outfilename = "gasvelocities.out"):
# 		super().__init__(onset = onset, dr = dr, dt = dt,
# 			outfilename = outfilename)
# 		self.maximum = maximum
# 		self.scale = scale

# 	def __call__(self, time):
# 		radii = [self.dr * i for i in range(int(MAX_RADIUS / self.dr))]
# 		vgas = [self.maximum * (1 - np.exp(-r / self.scale)) for r in radii]
# 		self.write(time, radii, vgas)
# 		return [radii, vgas]


class barflows:

	def __init__(self, maximum_speed,
		inneredge = 2, outeredge = 5, innerscale = 0.5, outerscale = 0.5):
		self.inner = logistic(minimum = 0, maximum = maximum_speed,
			midpoint = inneredge, scale = innerscale)
		self.outer = logistic(minimum = maximum_speed, maximum = 0,
			midpoint = outeredge, scale = outerscale)

	def __call__(self, r):
		return self.inner(r) * self.outer(r)


class oscillatory(base, sinusoid):

	def __init__(self, average, amplitude, period, phase = 0, onset = 1,
		dr = 0.1, dt = 0.01, outfilename = "gasvelocities.out"):
		base.__init__(self, onset = onset, dr = dr, dt = dt,
			outfilename = outfilename)
		sinusoid.__init__(self, amplitude = amplitude, period = period,
			phase = phase)
		self.average = average

	def __call__(self, time):
		radii = [self.dr * i for i in range(int(MAX_RADIUS / self.dr))]
		vgas = self.average + sinusoid.__call__(self, time)
		vgas = len(radii) * [vgas]
		self.write(time, radii, vgas)
		return [radii, vgas]



class linear(base):

	def __init__(self, dvdr = -0.1, onset = 1, dr = 0.1, dt = 0.01,
		outfilename = "gasvelocities.out"):
		super().__init__(onset = onset, dr = dr, dt = dt,
			outfilename = outfilename)
		self.dvdr = dvdr


	def __call__(self, time):
		radii = [self.dr * i for i in range(int(MAX_RADIUS / self.dr))]
		vgas = [self.dvdr * r for r in radii]
		# if self.bar is not None:
		# 	bareffect = [self.bar(r) for r in radii]
		# 	vgas = [a + b for a, b in zip(vgas, bareffect)]
		# else: pass
		self.write(time, radii, vgas)
		return [radii, vgas]



class potential_well_deepening(base):

	def __init__(self, mw_model, gamma = 0.2, onset = 1, dr = 0.1, dt = 0.01,
		recycling = 0.4, outfilename = "gasvelocities.out"):
		super().__init__(onset = onset, dr = dr, dt = dt,
			outfilename = outfilename)
		self.gamma = gamma
		self.mw_model = mw_model
		self.recycling = recycling
		self.evol = evoldata(mw_model, timestep = dt, recycling = recycling)


	def __call__(self, time):
		radii = [self.dr * i for i in range(self.mw_model.n_zones)]
		timestep = int(time / self.dt)
		sfr = 0
		mstar = 0
		for i in range(len(radii)):
			# decrement sfr by 1 - r because the potential well
			# deepening argument is based on the time derivative of the
			# stellar mass, which differs in important detail from the
			# specific star formation rate.
			sfr += (1 - self.recycling) * self.evol.sfrs[i][timestep]
			mstar += self.evol.mstars[i][timestep]
		vgas = [-r * self.gamma * sfr / mstar for r in radii]
		self.write(time, radii, vgas)
		return [radii, vgas]




class angular_momentum_dilution(base):


	def __init__(self, mw_model, beta_phi_in = 0.7, beta_phi_out = 0, onset = 1,
		dr = 0.1, dt = 0.01, outfilename = "gasvelocities.out"):
		super().__init__(onset = onset, dr = dr, dt = dt,
			outfilename = outfilename)
		self.mw_model = mw_model
		self.beta_phi_in = beta_phi_in
		self.beta_phi_out = beta_phi_out


	# def __call__(self, time, recycling = 0.4):
	# 	radii = [self.dr * i for i in range(int(MAX_RADIUS / self.dr))]
	# 	vgas = len(radii) * [0.]
	# 	crf = vice.cumulative_return_fraction(time)
	# 	def f(radius, vgas):
	# 		return self.dvdr(time, radius, vgas, recycling = crf)
	# 	profile = solve_ivp(f, [0, MAX_SF_RADIUS], [0], dense_output = True)
	# 	for i in range(1, len(radii)):
	# 		if radii[i] <= MAX_SF_RADIUS:
	# 			# vgas[i] = vgas[i - 1] + self.dr * self.dvdr(time, radii[i - 1],
	# 			# 	vgas[i - 1], recycling = recycling)
	# 			vgas[i] = profile.sol(radii[i])[0]
	# 			# print(vgas[i])
	# 		else:
	# 			vgas[i] = 0
	# 	self.write(time, radii, vgas)
	# 	return [radii, vgas]


	def __call__(self, time):
		radii = [self.dr * i for i in range(int(MAX_RADIUS / self.dr))]
		vgas = len(radii) * [0.]
		crf = vice.cumulative_return_fraction(time)
		for i in range(1, len(radii)):
			if radii[i] <= MAX_SF_RADIUS:
				vgas[i] = vgas[i - 1] + self.dr * self.dvdr(time, radii[i - 1],
					vgas[i - 1], recycling = crf)
			else:
				vgas[i] = 0
		self.write(time, radii, vgas)
		return [radii, vgas]


	def dvdr(self, time, radius, vgas, recycling = 0.4):
		zone = get_bin_number(self.mw_model.annuli, radius)
		if zone < 0: raise ValueError(
			"Radius outside of allowed range: %g" % (radius))

		neighbor = self.mw_model.zones[zone + 1]
		zone = self.mw_model.zones[zone]

		sfr = zone.func(time) * 1.e9 # yr^-1 -> Gyr^-1
		n_sfr = neighbor.func(time) * 1.e9
		sfr_next = zone.func(time + self.dt) * 1.e9

		taustar = zone.tau_star(time, sfr / 1.e9)
		n_taustar = neighbor.tau_star(time, n_sfr / 1.e9)
		taustar_next = zone.tau_star(time + self.dt, sfr_next / 1.e9)

		mg = sfr * taustar
		n_mg = n_sfr * n_taustar
		mg_next = sfr_next * taustar_next
		dlnmgdt = (mg_next - mg) / (mg * self.dt)

		sigmag = mg / (np.pi * ((radius + self.dr)**2 - radius**2))
		n_sigmag = n_mg / (np.pi * ((radius + 2 * self.dr)**2 -
			(radius + self.dr)**2))
		dlnsigmagdr = (n_sigmag - sigmag) / (sigmag * self.dr)

		if callable(zone.eta):
			eta = zone.eta(time)
		else:
			eta = zone.eta
		if callable(self.beta_phi_in):
			beta_phi_in = self.beta_phi_in(radius, time)
		else:
			beta_phi_in = self.beta_phi_in
		if callable(self.beta_phi_out):
			beta_phi_out = self.beta_phi_out(radius, time)
		else:
			beta_phi_out = self.beta_phi_out

		if radius:
			one_over_r = 1 / radius
		else:
			dlnmgdr = (n_mg - mg) / (mg * self.dr)
			one_over_r = dlnmgdr - dlnsigmagdr

		dvdr = 0
		dvdr -= dlnmgdt
		dvdr -= (1 - recycling) / taustar
		dvdr += eta / taustar * (beta_phi_out - beta_phi_in) / (beta_phi_in - 1)
		dvdr -= vgas * (one_over_r * (beta_phi_in - 2) / (beta_phi_in - 1)
			+ dlnsigmagdr)

		return dvdr




		# dvdr = 0

		# sfr = self.mw_model.zones[zone].func(time)
		# tau_star = self.mw_model.zones[zone].tau_star(time, sfr)
		# Mg = sfr * tau_star * 1.e9 # yr^-1 -> Gyr^-1
		# sfr_next = self.mw_model.zones[zone].func(time + self.dt)
		# Mg_next = sfr_next * self.mw_model.zones[zone].tau_star(
		# 	time + self.dt, sfr_next) * 1.e9
		# dlnMg_dt = (Mg_next - Mg) / (Mg * self.dt)
		# dvdr -= dlnMg_dt
		# dvdr -= (1 - recycling) / tau_star

		# if callable(self.mw_model.zones[zone].eta):
		# 	eta = self.mw_model.zones[zone].eta(time)
		# else:
		# 	eta = self.mw_model.zones[zone].eta
		# if callable(self.beta_phi_in):
		# 	beta_phi_in = self.beta_phi_in(radius, time)
		# else:
		# 	beta_phi_in = self.beta_phi_in
		# if callable(self.beta_phi_out):
		# 	beta_phi_out = self.beta_phi_out(radius, time)
		# else:
		# 	beta_phi_out = self.beta_phi_out

		# dvdr -= eta / tau_star * (1 - beta_phi_out) / (1 - beta_phi_in)

		# if radius + self.dr < MAX_SF_RADIUS:
		# 	Sigmag = Mg / (np.pi * (radius + self.dr)**2 - radius**2)
		# 	sfr_next = self.mw_model.zones[zone + 1].func(time)
		# 	Mg_next = sfr_next * self.mw_model.zones[zone + 1].tau_star(
		# 		time, sfr_next) * 1.e9 # yr^-1 -> Gyr^-1
		# 	Sigmag_next = Mg_next / (np.pi * (radius + 2 * self.dr)**2 - 
		# 		(radius + self.dr)**2)
		# 	dlnSigmag_dr = (Sigmag_next - Sigmag) / (Sigmag * self.dr)
		# else:
		# 	dlnSigmag_dr = -1 / self.dr # Sigmag_next -> 0

		# if radius:
		# 	dvdr -= vgas * (dlnSigmag_dr +
		# 		1 / radius * (beta_phi_in - 2) / (beta_phi_in - 1))
		# else:
		# 	# handle the 1 / r discontinuity by substituting in
		# 	# 1 / r = d\ln M / dr - d\ln\Sigma / dr -- algebraic solution below
		# 	dlnMg_dr = (Mg_next - Mg) / (Mg * self.dr)
		# 	dvdr -= vgas / (beta_phi_in - 1) * (
		# 		dlnSigmag_dr + (beta_phi_in - 2) * dlnMg_dr)

		# return dvdr




class river(base):

	def __init__(self, mw_model, onset = 1, dr = 0.1, dt = 0.01,
		outfilename = "gasvelocities.out"):
		super().__init__(onset = onset, dr = dr, dt = dt,
			outfilename = outfilename)
		if isinstance(mw_model, vice.milkyway):
			self.mw_model = mw_model
		else:
			raise TypeError(r"""\
Attribute 'mw_model' must be of type vice.milkyway. Got: %s.""" % (
				type(mw_model)))


	def __call__(self, time, recycling = 0.4):
		radii = [self.dr * i for i in range(int(MAX_RADIUS / self.dr))]
		vgas = len(radii) * [0.]
		for i in range(1, len(radii)):
			if radii[i] <= MAX_SF_RADIUS:
				crf = vice.cumulative_return_fraction(time)
				vgas[i] = self.next_vgas(vgas[i - 1], time, radii[i - 1],
					recycling = crf)
			else:
				vgas[i] = 0
		self.write(time, radii, vgas)
		return [radii, vgas]


	def next_vgas(self, vgas, time, radius, recycling = 0.4):
		zone = get_bin_number(self.mw_model.annuli, radius)
		if zone < 0: raise ValueError(
			"Radius outside of allowed range: %g" % (radius))

		radius = self.mw_model.annuli[zone] # use inner edge
		if zone == len(self.mw_model.annuli) - 1: return 0
		neighbor = self.mw_model.zones[zone + 1]
		zone = self.mw_model.zones[zone]

		sfr = zone.func(time) * 1.e9 # yr^-1 -> Gyr^-1
		n_sfr = neighbor.func(time) * 1.e9

		taustar = zone.tau_star(time, sfr / 1.e9)
		n_taustar = neighbor.tau_star(time, n_sfr / 1.e9)

		mgas = sfr * taustar
		n_mgas = n_sfr * n_taustar

		sfr_next = zone.func(time + self.dt) * 1.e9
		mgas_next = sfr_next * zone.tau_star(time + self.dt, sfr_next / 1.e9)

		if callable(zone.eta):
			eta = zone.eta(time)
		else:
			eta = zone.eta

		x = vgas**2 * self.dt**2 - 2 * radius * vgas * self.dt
		x /= 2 * radius * self.dr + self.dr**2
		x *= mgas
		x += mgas_next - mgas
		x += sfr * self.dt * (1 + eta - recycling) 
		x *= 2 * radius * self.dr + 3 * self.dr**2
		x /= n_mgas
		x += (radius + self.dr)**2

		if x < 0: raise ValueError("x < 0: %.5e. r = %.5e. t = %.5e" % (x,
			radius, time))

		n_vgas = 1 / self.dt * (radius + self.dr - np.sqrt(x))
		return n_vgas


# 	def __call__(self, time, recycling = 0.4, dr = 0.1, dt = 0.01):
# 		radii = [dr * i for i in range(int(MAX_RADIUS / dr))]
# 		vgas = len(radii) * [0.]
# 		vgas[0] = 0
# 		vgas[1] = self.v_at_deltaR(time, recycling = recycling, dr = dr, dt = dt)
# 		for i in range(2, len(radii)):
# 			if radii[i] <= MAX_SF_RADIUS:
# 				vgas[i] = vgas[i - 1] + dr * self.dvdr(time, radii[i - 1],
# 					vgas[i - 1], recycling = recycling, dr = dr, dt = dt)
# 			else:
# 				vgas[i] = 0
# 		self.write(time, radii, vgas)
# 		return [radii, vgas] # kpc / Gyr



# 	def v_at_deltaR(self, time, recycling = 0.4):
# 		if self.mw_model.mode != "sfr": raise ValueError("""\
# River model currently supports star formation mode.""")

# 		# SFR at R = [0, dR] at this timestep
# 		sfr = self.mw_model.zones[0].func(time) * 1.e9 # yr^-1 -> Gyr^-1
# 		# sfr = self.mw_model.zones[0].func(time)

# 		# SFR at R = [0, dR] at next timestep
# 		sfr_next = self.mw_model.zones[0].func(time + self.dt) * 1.e9
# 		# sfr_next = self.mw_model.zones[0].func(time + dt)

# 		# change in gas supply over this timestep
# 		dMg = sfr_next * self.mw_model.zones[0].tau_star(time + self.dt,
# 			sfr_next / 1e9)
# 		dMg -= sfr * self.mw_model.zones[0].tau_star(time, sfr / 1e9)
# 		# dMg = sfr_next * self.mw_model.zones[0].tau_star(time + dt, sfr_next)
# 		# dMg -= sfr * self.mw_model.zones[0].tau_star(time, sfr)
# 		# dMg *= 1.e9 # yr^-1 -> Gyr^-1 conversion in SFRs

# 		# gas supply at R = [dR, 2dR] at this timestep
# 		sfr_dr = self.mw_model.zones[1].func(time) * 1.e9
# 		Mg_dr = sfr_dr * self.mw_model.zones[1].tau_star(time, sfr_dr / 1e9)
# 		# Mg_dr *= 1.e9 # yr^-1 -> Gyr^-1 conversion in SFR

# 		if callable(self.mw_model.zones[0].eta):
# 			eta = self.mw_model.zones[0].eta(time)
# 		else:
# 			eta = self.mw_model.zones[0].eta
# 		x = (dMg + sfr * self.dt * (1 + eta - recycling)) / Mg_dr
# 		vgas = 1 - np.sqrt(1 + 3 * x)
# 		vgas *= self.dr / self.dt
# 		return vgas # kpc / Gyr


# 	def dvdr(self, time, radius, vgas, recycling = 0.4):
# 		zone = get_bin_number(self.mw_model.annuli, radius)
# 		if zone < 0: raise ValueError(
# 			"Radius outside of allowed range: %g" % (radius))
# 		else: pass

# 		# sfr = self.mw_model.zones[zone].func(time)
# 		# tau_star = self.mw_model.zones[zone].tau_star(time, sfr)
# 		sfr = self.mw_model.zones[zone].func(time) * 1.e9 # yr^-1 -> Gyr^-1
# 		tau_star = self.mw_model.zones[zone].tau_star(time, sfr / 1e9)
# 		# mgas = sfr * tau_star * 1.e9 # yr^-1 -> Gyr^-1 conversion in SFR
# 		mgas = sfr * tau_star
# 		# sfr_next = self.mw_model.zones[zone].func(time + dt)
# 		sfr_next = self.mw_model.zones[zone].func(time + self.dt) * 1e9
# 		# mgas_next = sfr_next * self.mw_model.zones[zone].tau_star(time + dt,
# 		# 	sfr_next) * 1.e9 # yr^-1 -> Gyr^-1 conversion in SFR
# 		mgas_next = sfr_next * self.mw_model.zones[zone].tau_star(time + self.dt,
# 			sfr_next / 1e9)
# 		dlnmgas_dt = (mgas_next - mgas) / (mgas * self.dt)

# 		# replace 1 / R + d\ln\Sigma_g / dR with d\ln M_g / dR
# 		if radius + self.dr < MAX_SF_RADIUS:
# 			# sfr_next = self.mw_model.zones[zone + 1].func(time)
# 			# mgas_next = sfr_next * self.mw_model.zones[zone + 1].tau_star(
# 			# 	time, sfr_next) * 1.e9 # yr^-1 -> Gyr^-1 conversion in SFR
# 			sfr_next = self.mw_model.zones[zone + 1].func(time) * 1.e9
# 			mgas_next = sfr_next * self.mw_model.zones[zone + 1].tau_star(
# 				time, sfr_next / 1e9)
# 			dlnmgas_dr = (mgas_next - mgas) / (mgas * self.dr)
# 		else:
# 			# this is the last zone and the calculation is stopping here anyway
# 			dlnmgas_dr = -1 / self.dr # mgas_next -> 0

# 		if callable(self.mw_model.zones[0].eta):
# 			eta = self.mw_model.zones[0].eta(time)
# 		else:
# 			eta = self.mw_model.zones[0].eta

# 		dvdr = -dlnmgas_dt
# 		dvdr -= (1 + eta - recycling) / tau_star
# 		dvdr -= vgas * dlnmgas_dr
# 		return dvdr # Gyr^-1


