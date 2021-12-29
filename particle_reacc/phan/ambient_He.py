import numpy as np
import astropy.units as un
from astropy.constants import c, m_p, m_e, u

def ambient_He(ene_He):

    '''
    Cosmic ray spectrum for He
    :param ene: kinetic energy list (astropy quantity in GeV)
    '''
    # Check if input kinetic energy is GeV
    if not ene_He.unit is un.GeV:
        ene_He = ene_He.to(un.GeV)
    
    cc = c.to(un.cm/un.s)
    ene = ene_He.value
    PA = 4.002
    m_He = PA * u
    E_He = m_He * c**2
    rest_He = E_He.to('GeV')
    rest = rest_He.value

    def Beta(ene):
        gamma = (ene + rest)/rest
        Beta = np.sqrt(1-gamma**(-2))
        return Beta

    # Parameters
    const = 1*4*np.pi*195.4/un.m**2/un.GeV/un.s
    a = 1.03
    b = 3.18
    d = 1.21
    k = 0.77
    # He distribution
    
    n_p = (const/c) * (ene**a)/Beta(ene)**3 * ((ene**d + k**d)/(1+k**d))**(-b)
    
    return n_p.to(un.eV**(-1)/un.cm**3)

