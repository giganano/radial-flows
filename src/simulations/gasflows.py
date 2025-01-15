r"""
Handles radial gas flows in these models.
"""

from .._globals import MAX_SF_RADIUS, END_TIME
from .models.utils import get_bin_number
from vice.toolkit.interpolation import interp_scheme_1d
from vice.milkyway.milkyway import _MAX_RADIUS_ as MAX_RADIUS # 20 kpc
import numpy as np
import warnings
import vice


class driver(interp_scheme_1d):

	def __init__(self, *args, dt = 0.01, **kwargs):
		super().__init__(*args, **kwargs)
		self.dt = dt
	
	def __call__(self, x):
		test = super().__call__(x) * self.dt / 0.01
		if test < 0:
			return 0
		elif test > 1:
			return 1
		else:
			return test



class base:

	NORMTIME = 0.01 # Gyr

	def __init__(self, onset = 1, outfilename = "gasvelocities.out"):
		if outfilename is not None:
			self.outfile = open(outfilename, "w")
			self.outfile.write("# Time [Gyr]    ")
			self.outfile.write("Radius [kpc]    ")
			self.outfile.write("ISM radial velocity [kpc/Gyr]\n")
		else:
			self.outfile = None
		self.onset = onset


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


	def setup(self, mw_model, dr = 0.1, dt = 0.01, **kwargs):
		vgas_alltimes = []
		n_zones = int(MAX_RADIUS / dr)
		times = [dt * i for i in range(int(END_TIME / dt) + 10)]
		for i in range(len(times)):
			if i > self.onset / dt:
				radii, vgas = self.__call__(i * dt, **kwargs)
				vgas_alltimes.append(vgas)
			else:
				radii = [dr * i for i in range(n_zones)]
				vgas = len(radii) * [0.]
				vgas_alltimes.append(vgas)
		matrix_elements_inward = []
		matrix_elements_outward = []
		for i in range(n_zones):
			areafracs_inward = []
			areafracs_outward = []
			vgas = [row[i] for row in vgas_alltimes]
			for j in range(len(times)):
				if j > self.onset / dt:
					radius = i * dr
					if vgas[j] > 0: # outward flow
						numerator = 2 * (radius + dr) * vgas[j] * self.NORMTIME
						numerator -= vgas[j]**2 * self.NORMTIME**2
					else: # inward flow
						numerator = vgas[j]**2 * 0.01**2
						numerator -= 2 * radius * vgas[j] * self.NORMTIME
					denominator = 2 * radius * dr + dr**2
					areafrac = numerator / denominator
					if areafrac * dt / self.NORMTIME > 1:
						warnings.warn("""\
Area fraction larger than 1. Consider comparing results with different \
timestep sizes to assess the impact of numerical artifacts.""")
						areafrac = self.NORMTIME / dt - 1.e-9
					elif areafrac < 0:
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
				driver(times, areafracs_outward, dt = dt))
			matrix_elements_inward.append(
				driver(times, areafracs_inward, dt = dt))
		for i in range(n_zones):
			for j in range(n_zones):
				if i - 1 == j: # inward flows
					mw_model.migration.gas[i][j] = matrix_elements_inward[i]
				elif i + 1 == j: # outward flows
					mw_model.migration.gas[i][j] = matrix_elements_outward[i]
				else:
					mw_model.migration.gas[i][j] = 0




class constant(base):

	def __init__(self, speed, onset = 1, outfilename = "gasvelocities.out"):
		super().__init__(onset = onset, outfilename = outfilename)
		self.speed = speed


	def __call__(self, time, dr = 0.1, dt = 0.01, **kwargs):
		radii = [dr * i for i in range(int(MAX_RADIUS / dr))]
		if callable(self.speed):
			# it's a constant in radius, but not necessarily in time
			speed = self.speed(time, **kwargs)
		else:
			speed = self.speed
		vgas = len(radii) * [speed]
		self.write(time, radii, vgas)
		return [radii, vgas]



class linear(base):

	def __init__(self, dvdr = -0.1, onset = 1,
		outfilename = "gasvelocities.out"):
		super().__init__(onset = onset, outfilename = outfilename)
		self.dvdr = dvdr


	def __call__(self, time, dr = 0.1, dt = 0.01):
		radii = [dr * i for i in range(int(MAX_RADIUS / dr))]
		vgas = [self.dvdr * r for r in radii]
		self.write(time, radii, vgas)
		return [radii, vgas]




class angular_momentum_dilution(base):


	def __init__(self, mw_model, beta_phi_in = 0.7, beta_phi_out = 0, onset = 1,
		outfilename = "gasvelocities.out"):
		super().__init__(onset = onset, outfilename = outfilename)
		self.mw_model = mw_model
		self.beta_phi_in = beta_phi_in
		self.beta_phi_out = beta_phi_out


	def __call__(self, time, recycling = 0.4, dr = 0.1, dt = 0.01):
		radii = [dr * i for i in range(int(MAX_RADIUS / dr))]
		vgas = len(radii) * [0.]
		for i in range(1, len(radii)):
			if radii[i] <= MAX_SF_RADIUS:
				vgas[i] = vgas[i - 1] + dr * self.dvdr(time, radii[i - 1],
					vgas[i - 1], recycling = recycling, dr = dr, dt = dt)
			else:
				vgas[i] = 0
		self.write(time, radii, vgas)
		return [radii, vgas]


	def dvdr(self, time, radius, vgas, recycling = 0.4, dr = 0.1, dt = 0.01):
		zone = get_bin_number(self.mw_model.annuli, radius)
		if zone < 0: raise ValueError(
			"Radius outside of allowed range: %g" % (radius))
		dvdr = 0

		sfr = self.mw_model.zones[zone].func(time)
		tau_star = self.mw_model.zones[zone].tau_star(time, sfr)
		Mg = sfr * tau_star * 1.e9 # yr^-1 -> Gyr^-1
		sfr_next = self.mw_model.zones[zone].func(time + dt)
		Mg_next = sfr_next * self.mw_model.zones[zone].tau_star(
			time + dt, sfr_next) * 1.e9
		dlnMg_dt = (Mg_next - Mg) / (Mg * dt)
		dvdr -= dlnMg_dt
		dvdr -= (1 - recycling) / tau_star

		if callable(self.mw_model.zones[zone].eta):
			eta = self.mw_model.zones[zone].eta(time)
		else:
			eta = self.mw_model.zones[zone].eta
		if callable(self.beta_phi_in):
			beta_phi_in = self.beta_phi_in(radius, time)
		else:
			beta_phi_in = self.beta_phi_in
		if callable(self.beta_phi_out):
			beta_phi_out = self.beta_phi_out(radius, time)
		else:
			beta_phi_out = self.beta_phi_out

		dvdr -= eta / tau_star * (1 - beta_phi_out) / (1 - beta_phi_in)

		if radius + dr < MAX_SF_RADIUS:
			Sigmag = Mg / (np.pi * (radius + dr)**2 - radius**2)
			sfr_next = self.mw_model.zones[zone + 1].func(time)
			Mg_next = sfr_next * self.mw_model.zones[zone + 1].tau_star(
				time, sfr_next) * 1.e9 # yr^-1 -> Gyr^-1
			Sigmag_next = Mg_next / (np.pi * (radius + 2 * dr)**2 - 
				(radius + dr)**2)
			dlnSigmag_dr = (Sigmag_next - Sigmag) / (Sigmag * dr)
		else:
			dlnSigmag_dr = -1 / dr # Sigmag_next -> 0

		if radius:
			dvdr -= vgas * (dlnSigmag_dr +
				1 / radius * (beta_phi_in - 2) / (beta_phi_in - 1))
		else:
			# handle the 1 / r discontinuity by substituting in
			# 1 / r = d\ln M / dr - d\ln\Sigma / dr -- algebraic solution below
			dlnMg_dr = (Mg_next - Mg) / (Mg * dr)
			dvdr -= vgas / (beta_phi_in - 1) * (
				dlnSigmag_dr + (beta_phi_in - 2) * dlnMg_dr)

		return dvdr





class river(base):

	# EXTRA_RESOLUTION = 10

	def __init__(self, mw_model, onset = 1, outfilename = "gasvelocities.out"):
		super().__init__(onset = onset, outfilename = outfilename)
		if isinstance(mw_model, vice.milkyway):
			self.mw_model = mw_model
		else:
			raise TypeError(r"""\
Attribute 'mw_model' must be of type vice.milkyway. Got: %s.""" % (
				type(mw_model)))


	def __call__(self, time, recycling = 0.4, dr = 0.1, dt = 0.01):
		radii = [dr * i for i in range(int(MAX_RADIUS / dr))]
		vgas = len(radii) * [0.]
		vgas[0] = 0
		vgas[1] = self.v_at_deltaR(time, recycling = recycling, dr = dr, dt = dt)
		for i in range(2, len(radii)):
			if radii[i] <= MAX_SF_RADIUS:
				vgas[i] = vgas[i - 1] + dr * self.dvdr(time, radii[i - 1],
					vgas[i - 1], recycling = recycling, dr = dr, dt = dt)
			else:
				vgas[i] = 0
		self.write(time, radii, vgas)
		return [radii, vgas] # kpc / Gyr


	def v_at_deltaR(self, time, recycling = 0.4, dr = 0.1, dt = 0.01):
		if self.mw_model.mode != "sfr": raise ValueError("""\
River model currently supports star formation mode.""")

		# SFR at R = [0, dR] at this timestep
		sfr = self.mw_model.zones[0].func(time)

		# SFR at R = [0, dR] at next timestep
		sfr_next = self.mw_model.zones[0].func(time + dt)

		# change in gas supply over this timestep
		dMg = sfr_next * self.mw_model.zones[0].tau_star(time + dt, sfr_next)
		dMg -= sfr * self.mw_model.zones[0].tau_star(time, sfr)
		dMg *= 1.e9 # yr^-1 -> Gyr^-1 conversion in SFRs

		# gas supply at R = [dR, 2dR] at this timestep
		sfr_dr = self.mw_model.zones[1].func(time)
		Mg_dr = sfr_dr * self.mw_model.zones[1].tau_star(time, sfr_dr)
		Mg_dr *= 1.e9 # yr^-1 -> Gyr^-1 conversion in SFR

		if callable(self.mw_model.zones[0].eta):
			eta = self.mw_model.zones[0].eta(time)
		else:
			eta = self.mw_model.zones[0].eta
		x = (dMg + 1.e9 * sfr * dt * (1 + eta - recycling)) / Mg_dr
		vgas = 1 - np.sqrt(1 + 3 * x)
		vgas *= dr / dt
		# if abs(time - 0.11) < 1.e-3:
		# 	print("============================================")
		# 	print(dMg + 1.e9 * sfr * dt * (1 + eta - recycling))
		# 	print(Mg_dr)
		# 	print(vgas)
		# 	print(Mg_dr * (vgas**2 * dt**2 - 2 * dr * vgas * dt) / (3 * dr**2))
		# 	print("============================================")
		# else: pass
		return vgas # kpc / Gyr


	def dvdr(self, time, radius, vgas, recycling = 0.4, dr = 0.1, dt = 0.01):
		zone = get_bin_number(self.mw_model.annuli, radius)
		if zone < 0: raise ValueError(
			"Radius outside of allowed range: %g" % (radius))
		else: pass

		sfr = self.mw_model.zones[zone].func(time)
		tau_star = self.mw_model.zones[zone].tau_star(time, sfr)
		mgas = sfr * tau_star * 1.e9 # yr^-1 -> Gyr^-1 conversion in SFR
		sfr_next = self.mw_model.zones[zone].func(time + dt)
		mgas_next = sfr_next * self.mw_model.zones[zone].tau_star(time + dt,
			sfr_next) * 1.e9 # yr^-1 -> Gyr^-1 conversion in SFR
		dlnmgas_dt = (mgas_next - mgas) / (mgas * dt)

		# replace 1 / R + d\ln\Sigma_g / dR with d\ln M_g / dR
		if radius + dr < MAX_SF_RADIUS:
			sfr_next = self.mw_model.zones[zone + 1].func(time)
			mgas_next = sfr_next * self.mw_model.zones[zone + 1].tau_star(
				time, sfr_next) * 1.e9 # yr^-1 -> Gyr^-1 conversion in SFR
			dlnmgas_dr = (mgas_next - mgas) / (mgas * dr)
		else:
			# this is the last zone and the calculation is stopping here anyway
			dlnmgas_dr = -1 / dr # mgas_next -> 0

		if callable(self.mw_model.zones[0].eta):
			eta = self.mw_model.zones[0].eta(time)
		else:
			eta = self.mw_model.zones[0].eta

		dvdr = -dlnmgas_dt
		dvdr -= (1 + eta - recycling) / tau_star
		dvdr -= vgas * dlnmgas_dr
		return dvdr # Gyr^-1


