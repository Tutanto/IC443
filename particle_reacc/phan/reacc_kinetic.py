import numpy as np
import astropy.units as u
from astropy.constants import c, m_p, m_e
from scipy import integrate
from ambient_CR_kinetic import ambient_CR_kinetic
from distr2spectr_kinetic import distr2spectr_kinetic

def reacc_kinetic(ene,co,br,s,particle):

    # Check if input energy is MeV
    if not ene.unit is u.MeV or not co.unit is u.MeV or not br.unit is u.MeV:
        ene = ene.to(u.MeV)
        co = co.to(u.MeV)
        br = br.to(u.MeV)
#        raise Exception("Please use MeV unit")

    r = 4
    alpha = (r + 2)/(r - 1) 
    ene = ene*s**(-1/3)
    ene_1 = ene[ene < br].value
    ene_2 = ene[ene > br].value
    cc = c.to(u.cm/u.s)
    co = co.value
    br = br.value
    if particle == 'proton':
        const = 4*np.pi*1.882e-9*u.s/u.cm**4/u.eV
        i1 = 0.129
        i2 = 2.829
        E_br = 624.5
        E_p = m_p * c**2
        rest = E_p.to('MeV')
        rest = rest.value

    if particle == 'electron':
        const = 4*np.pi*4.658e-7*u.s/u.cm**4/u.eV
        i1 = -1.236
        i2 = 2.033
        E_br = 736.2
        E_e = m_e * c**2
        rest = E_e.to('MeV')
        rest = rest.value
        
    def beta(ene):
        gamma = ((ene+rest)/rest)
        beta = np.sqrt(1-(1/gamma)**2)
        return beta

    co = np.sqrt(co*(co+2*rest))/cc
    br = np.sqrt(br*(br+2*rest))/cc

    f = lambda x: ((np.sqrt(x*(x+2*rest))/cc.value)**(alpha-1)) * (x**i1 * (1+x/E_br)**(-i2))/(beta(x)*cc.value)


    np.vectorize(f)

    def fint(a,z):
        return integrate.quad(f, a, z)

    fint=np.vectorize(fint)

    an_1=fint(min(ene.value), ene_1)
    an_2=fint(min(ene.value), ene_2)

    def p2k(ene):
        p = (np.sqrt(ene**2 + 2*rest*ene))/cc
        return p
    def dp_dk(ene):
        var = (ene + rest)/(cc*np.sqrt(ene*(ene+2*rest)))
        return var
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    n_acc_1 = an_1[0] * (alpha+2) * (p2k(ene_1))**(-alpha) * np.exp(-(p2k(ene_1))/co) * dp_dk(ene_1) * const
    n_acc_2 = an_2[0] * (alpha+2) * (p2k(ene_2))**(-alpha) * br/(p2k(ene_2)) * np.exp(-(p2k(ene_2))/co) * dp_dk(ene_2) * const

    f_p, f_e = ambient_CR_kinetic(ene,ene)
    n_p, n_e = distr2spectr_kinetic(f_p,f_e,ene,ene)
    
    # Check if n_acc as same units as n_p
    if not n_acc_1.unit == n_p.unit or not n_acc_2.unit == n_p.unit:
        print(n_acc_1.unit)
        print(n_p.unit)
        raise Exception("***WARNING*** n_acc and n_p must have same units")


    idx = find_nearest(ene,ene_2[0])

    if particle == 'proton':
        n_fin = s**(2/3) * np.concatenate((n_acc_1.value, np.maximum(n_p.value[idx:],n_acc_2.value)))
    if particle == 'electron':
        n_fin = s**(2/3) * np.concatenate((n_acc_1.value, np.maximum(n_e.value[idx:],n_acc_2.value)))
    return n_fin
