import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.constants import c, m_p, m_e
from reacc_kinetic import reacc_kinetic
from ambient_CR_kinetic import ambient_CR_kinetic
from distr2spectr_kinetic import distr2spectr_kinetic

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def Beta_p(ene):
    E_p = m_p * c**2
    rest = E_p.to('MeV')
    rest = rest.value
    gamma = ((ene+rest)/rest)
    Beta = np.sqrt(1-(1/gamma)**2)
    return Beta

def Beta_e(ene):
    E_e = m_e * c**2
    rest = E_e.to('MeV')
    rest = rest.value
    gamma = ((ene+rest)/rest)
    Beta = np.sqrt(1-(1/gamma)**2)
    return Beta

# Prepare an energy array for saving the particle distribution
energy = np.logspace(-3, 7, 100) * u.MeV

cc = c.to(u.cm/u.s)
rest_p = (m_p * c**2).to(u.MeV)
rest_e = (m_e * c**2).to(u.MeV)

co = 100 * u.GeV
br = 50 * u.GeV
s = 12

j_p,j_e = ambient_CR_kinetic(12**(-1/3)*energy,12**(-1/3)*energy)
n_p,n_e = distr2spectr_kinetic(j_p,j_e,12**(-1/3)*energy,12**(-1/3)*energy)

n_pr_reacc = reacc_kinetic(energy,co.to(u.MeV),br.to(u.MeV),s,particle='proton')
n_ele_reacc = reacc_kinetic(energy,co.to(u.MeV),br.to(u.MeV),s,particle='electron')

idx = find_nearest(energy,1.e3)
#n_pr_reacc = n_pr_reacc/n_pr_reacc[idx]
#n_ele_reacc = n_ele_reacc/n_ele_reacc[idx]


#plt.loglog(energy,j_p,'k',label='j$_{p} (Vojager)$')
#plt.loglog(energy,j_e,'r',label='j$_{e} (Vojager)$')
plt.loglog(energy,s**(2/3)*n_p*(energy**2),'k--',label='proton')
plt.loglog(energy,s**(2/3)*n_e*(energy**2),'r--',label='electron')
'''plt.xlabel('K [MeV]')
#plt.ylabel('j(K) [cm$^{-2}$eV$^{-1}$sr$^{-1}$s$^{-1}$]')
#plt.ylabel('dN/dK')
plt.ylabel('K$^2$ dN/dK')
plt.legend()
plt.show()
'''
plt.loglog(energy,n_pr_reacc*((energy)**2),'k',label='n$_{ap} (proton)$')
plt.loglog(energy,n_ele_reacc*((energy)**2),'r',label='n$_{ap} (electron)$')
plt.xlabel('K [MeV]')
plt.ylabel('K$^2$ dN/dK')
#plt.ylim(1.e-2,1.e2)
#plt.xlim(1.e0,1.e6)
plt.legend()
plt.show()

