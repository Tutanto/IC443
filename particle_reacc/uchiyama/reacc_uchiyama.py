import numpy as np
import astropy.units as u
from astropy.constants import c, m_p, m_e
from scipy import integrate
from ambient_CR_uchiyama import ambient_CR_uchiyama

def reacc_uchiyama(ene,co,br,s,particle):

    # Check if input energy is GeV
    if not ene.unit is u.GeV or not co.unit is u.GeV or not br.unit is u.GeV:
        ene = ene.to(u.GeV)
        co = co.to(u.GeV)
        br = br.to(u.GeV)

    r = 4
    alpha = (r + 2)/(r - 1) 
    ene = ene*s**(-1/3)
    ene_1 = ene[ene < br].value
    ene_2 = ene[ene > br].value
    cc = c.to(u.cm/u.s)
    co = co.value
    br = br.value
    if particle == 'proton':
        E_p = m_p * c**2
        rest = E_p.to('GeV')
        rest = rest.value
        const = 4*np.pi*1.9*u.s/u.cm**4/u.GeV
        a = 1.5
        b = -2.76
        def beta(mom):
            factor = mom/rest
            beta = factor/np.sqrt(1+(factor)**2)
            return beta


        f = lambda x: ((np.sqrt(x*(x+2*rest))/cc.value)**(alpha-1)) * beta(np.sqrt(x**2 + 2*rest*x))**a * (np.sqrt(x**2 + 2*rest*x))**b * (x + rest)/(np.sqrt(x*(x+2*rest))*cc.value)

    if particle == 'electron':
        E_e = m_e * c**2
        rest = E_e.to('GeV')
        rest = rest.value
        const = 4*np.pi*2e-2*u.s/u.cm**4/u.GeV
        a = 2
        b = -0.55

        f = lambda x: ((np.sqrt(x*(x+2*rest))/cc.value)**(alpha-1)) * (np.sqrt(x**2 + 2*rest*x))**(-a) * (1 + (np.sqrt(x**2 + 2*rest*x))**a)**b * (x + rest)/(np.sqrt(x*(x+2*rest))*cc.value)

        
    co = np.sqrt(co*(co+2*rest))/cc
    br = np.sqrt(br*(br+2*rest))/cc

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
    
    n_acc_1 = n_acc_1.to(u.eV**(-1)/u.cm**3)
    n_acc_2 = n_acc_2.to(u.eV**(-1)/u.cm**3)

    n_p, n_e = ambient_CR_uchiyama(ene,ene)
    n_p = n_p.to(u.eV**(-1)/u.cm**3)
    n_e = n_e.to(u.eV**(-1)/u.cm**3)

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
