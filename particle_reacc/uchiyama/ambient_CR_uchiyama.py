import numpy as np
import astropy.units as u
from astropy.constants import c, m_p, m_e

def ambient_CR_uchiyama(ene_p,ene_e):

    '''
    Cosmic ray spectrum for proton ad electrons
    Uchiyama et al. (2010)
    :param ene: proton energy list (astropy quantity in GeV)
    '''
    # Check if input kinetic energy is GeV
    if not ene_p.unit is u.GeV:
        ene_p = ene_p.to(u.GeV)
    
    cc = c.to(u.cm/u.s)
    ene = ene_p.value
    E_p = m_p * c**2
    rest = E_p.to('GeV')
    rest_p = rest.value

    def Beta_p(mom):
        factor = mom*cc.value/rest_p
        Beta = factor/np.sqrt(1+factor**2) 
        return Beta

    # Parameters protons
    const = 4*np.pi*1.9/u.cm**2/u.GeV/u.s
    a = 1.5
    b = -2.76
    # proton distribution
    mom = np.sqrt(ene**2 + 2*rest_p*ene)#/cc.value
    dp_dk = (ene + rest_p)/(cc*np.sqrt(ene*(ene+2*rest_p)))
    
    n_p = const * Beta_p(mom)**a * mom**b * dp_dk

    # Check if input kinetic energy is GeV
    if not ene_e.unit is u.GeV:
        ene_e = ene_e.to(u.GeV)
 
    ene = ene_e.value
    E_e = m_e * c**2
    rest = E_e.to('GeV')
    rest_e = rest.value

    # Parameters electrons
    const = 4*np.pi*2.e-2/u.cm**2/u.GeV/u.s
    a = 2
    b = -0.55
    # electron distribution
    mom = np.sqrt(ene**2 + 2*rest_e*ene)#/cc.value
    dp_dk = (ene + rest_e)/(cc*np.sqrt(ene*(ene+2*rest_e)))

    n_e = const*mom**(-a) * (1+mom**a)**b * dp_dk

    return n_p, n_e

