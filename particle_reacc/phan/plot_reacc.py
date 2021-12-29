import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.constants import c, m_p, m_e
from reacc_kinetic import reacc_kinetic
from reacc_He import reacc_He

# Prepare an energy array for saving the particle distribution
energy = np.logspace(-3, 7, 100) * u.MeV

co = 100 * u.GeV
br = 50 * u.GeV
s = 12

n_pr = reacc_kinetic(energy,co,br,s,particle='proton') /u.cm**3/u.eV
n_he = reacc_He(energy,co,br,s) /u.cm**3/u.eV


plt.loglog(energy.to(u.MeV),n_pr*energy.to(u.eV)**2,'k',label='n$_{ap}$ (proton Phan)')
plt.loglog(energy.to(u.MeV),n_he*energy.to(u.eV)**2,'r',label='n$_{ap}$ (He)')
plt.xlabel('K [MeV/n]')
plt.ylabel('K$^2$ dN/dK [eV cm$^{-3}$]')
plt.ylim(1.e-4,1.e3)
plt.xlim(1.e0,1.e7)
plt.legend()
plt.show()

