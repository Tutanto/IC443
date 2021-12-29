import numpy as np
import astropy.units as u
from astropy.constants import c, m_p, m_e

def ambient_CR(ene_p):

    '''
    Cosmic ray spectrum for protons in kinetic energy
    :param ene: proton energy list (astropy quantity in GeV)
    '''
    # Check if input kinetic energy is GeV
    if not ene_p.unit is u.GeV:
        ene_p = ene_p.to(u.GeV)
    
    ene = ene_p.value
    E_p = m_p * c**2
    rest = E_p.to('GeV')
    rest_p = rest.value

    def Beta_p(ene):
        factor = (ene+rest_p)/rest_p
        Beta = factor/np.sqrt(1+factor**2) 
        return Beta

    # Parameters protons
    const = 1.56e-10/u.cm**3/u.GeV
    a = 1.03
    b = 1.21
    d = -3.18
    # proton distribution
    n_p = const * (ene**a)/Beta_p(ene)**3 * ((ene**a + 0.77**b)/(1 + 0.77**b))**d

    return n_p

