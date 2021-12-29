import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.constants import c, m_p, m_e
from reacc import reacc

# Prepare an energy array for saving the particle distribution
energy = np.logspace(-3, 7, 100) * u.MeV

co = 100 * u.GeV
br = 50 * u.GeV
s = 12

n_pr_reacc = reacc(energy,co,br,s) /u.cm**3/u.eV


plt.loglog(energy,n_pr_reacc*((energy.to(u.eV))**2),'k',label='n$_{ap} (proton)$')
plt.xlabel('K [MeV]')
plt.ylabel('K$^2$ dN/dK [eV cm$^{-3}$]')
plt.ylim(1.e-4,1.e5)
plt.xlim(1.e0,1.e7)
plt.legend()
plt.show()

