import numpy as np
import astropy.units as u
from astropy.constants import c, m_p

def ambient_CR_kinetic(ene_p,ene_e):

    '''
    Cosmic ray intensity for proton ad electrons
    V. H. M. Phan et al. (2018)
    :param ene: proton energy list (astropy quantity in MeV)
    '''
    # Check if input kinetic energy is MeV
    if not ene_p.unit is u.MeV:
        raise Exception("Please use MeV unit for energy")
    
    ene = ene_p.value

    # Parameters protons
    const = 1.882e-9/u.cm**2/u.eV/u.s/u.sr
    a = 0.129
    b = 2.829
    E_br = 624.5
    # proton distribution
    j_p = const*ene**a * (1+ene/E_br)**(-b)

    # Check if input energy is MeV
    if not ene_e.unit is u.MeV:
        raise Exception("Please use MeV unit for energy")

    ene = ene_e.value

    # Parameters electrons
    const = 4.658e-7/u.cm**2/u.eV/u.s/u.sr
    a = -1.236
    b = 2.033
    E_br = 736.2
    # electron distribution
    j_e = const*ene**a * (1+ene/E_br)**(-b)

    return j_p, j_e

