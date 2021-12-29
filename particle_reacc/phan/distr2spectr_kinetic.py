import numpy as np
import astropy.units as u
from astropy.constants import c, m_p, m_e

def distr2spectr_kinetic(j_p,j_e,ene_p,ene_e):

    '''
    Create number density distribution
    of ambient Cosmic Ray proton and electron
    assuming particle distribution from ambient_CR
    :param ene: proton energy list (astropy quantity in MeV)
    :param f: particle distribution
    '''
    # Check if input energy is MeV
    if not ene_p.unit is u.MeV:
        raise Exception("Please use MeV unit for energy")

    ene_p = ene_p.value
    # Compute beta
    def Beta_p(ene):
        E_p = m_p * c**2
        rest = E_p.to('MeV')
        rest = rest.value
        gamma = ((ene+rest)/rest)
        Beta = np.sqrt(1-(1/gamma)**2)
        return Beta
    n_p = 4*np.pi*u.sr*j_p/(Beta_p(ene_p)*c)
    n_p = n_p.to(u.eV**(-1)*u.cm**(-3))
    
    # Check if input energy is MeV
    if not ene_e.unit is u.MeV:
        raise Exception("Please use MeV unit for energy")

    ene_e = ene_e.value
    # Compute beta
    def Beta_e(ene):
        E_e = m_e * c**2
        rest = E_e.to('MeV')
        rest = rest.value
        gamma = ((ene+rest)/rest)
        Beta = np.sqrt(1-(1/gamma)**2)
        return Beta
    n_e = 4*np.pi*u.sr*j_e/(Beta_e(ene_e)*c)
    n_e = n_e.to(u.eV**(-1)*u.cm**(-3))


    return n_p, n_e

