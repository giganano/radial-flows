r"""
Handles radial gas flows in these models.
"""

from .._globals import MAX_SF_RADIUS
from .models.utils import get_bin_number
from vice.milkyway.milkyway import _MAX_RADIUS_ as MAX_RADIUS # 20 kpc
import numpy as np
import vice


class base:

	def __init__(self, outfilename = "gasvelocities.out"):
		if outfilename is not None:
			self.outfile = open(outfile, "w")
		else:
			self.outfile = None


	def __enter__(self):
		return self


	def __exit__(self, exc_type, exc_value, exc_tb):
		self.outfile.close()
		return exc_value is None


	def write(self, time, radii, velocities):
		if self.outfile is not None:
			for i in range(len(radii)):
				self.outfile.write("%.2e\t%.2e\t%.2e\n" % (
					time, radii[i], velocities[i]))
		else: pass


class river(base):

	def __init__(self, mw_model, outfilename = "gasvelocities.out"):
		super().__init__(outfilename = outfilename)
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

		if callable(self.mw_model[0].eta):
			eta = self.mw_model[0].eta(time)
		else:
			eta = self.mw_model[0].eta
		x = (dMg + 1.e9 * sfr * dt * (1 + eta - recycling)) / Mg_dr
		vgas = 1 - np.sqrt(1 + 3 * x)
		vgas *= dr / dt
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
			dlnmgas_dr = -1 / dr # (mgas_next -> 0)

		if callable(self.mw_model[0].eta):
			eta = self.mw_model[0].eta(time)
		else:
			eta = self.mw_model[0].eta

		dvdr = -dlnmgas_dt
		dvdr -= (1 + eta - recycling) / tau_star
		dvdr -= vgas * dlnmgas_dr
		return dvdr # Gyr^-1


