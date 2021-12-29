import numpy as np
import astropy.units as un
from astropy.constants import c, m_p, m_e, u
from scipy import integrate
from ambient_He import ambient_He

def reacc_He(ene,co,br,s):

    # Check if input energy is GeV
    if not ene.unit is un.GeV or not co.unit is un.GeV or not br.unit is un.GeV:
        ene = ene.to(un.GeV)
        co = co.to(un.GeV)
        br = br.to(un.GeV)

    r = 4
    alpha = (r + 2)/(r - 1) 
    ene = ene*s**(-1/3)
    ene_1 = ene[ene < br].value
    ene_2 = ene[ene > br].value
    cc = c.to(un.cm/un.s)
    co = co.value
    br = br.value
    
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


    f = lambda x: ((np.sqrt(x*(x+2*rest))/cc.value)**(alpha-1)) * (x**a)/(cc.value*Beta(x)**3) * ((x**d + k**d)/(1 + k**d))**(-b) * (x + rest)/(np.sqrt(x*(x+2*rest)))

    np.vectorize(f)

    def fint(a,z):
        return integrate.quad(f, a, z)

    fint=np.vectorize(fint)

    an_1=fint(min(ene.value), ene_1)
    an_2=fint(min(ene.value), ene_2)

    def p2k(ene):
        p = (np.sqrt(ene**2 + 2*rest*ene))/cc.value
        return p
    def dp_dk(ene):
        var = (ene + rest)/(cc*np.sqrt(ene*(ene+2*rest)))
        return var
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx
    
    co = np.sqrt(co*(co+2*rest))/cc.value
    br = np.sqrt(br*(br+2*rest))/cc.value

    n_acc_1 = an_1[0] * (alpha+2) * (p2k(ene_1))**(-alpha) * np.exp(-(p2k(ene_1))/co) * dp_dk(ene_1) * const
    n_acc_2 = an_2[0] * (alpha+2) * (p2k(ene_2))**(-alpha) * br/(p2k(ene_2)) * np.exp(-(p2k(ene_2))/co) * dp_dk(ene_2) * const 
    
    n_acc_1 = n_acc_1.to(un.eV**(-1)/un.cm**3)
    n_acc_2 = n_acc_2.to(un.eV**(-1)/un.cm**3)

    n_p  = ambient_He(ene)
    n_p = n_p.to(un.eV**(-1)/un.cm**3)

    # Check if n_acc as same units as n_p
    if not n_acc_1.unit == n_p.unit or not n_acc_2.unit == n_p.unit:
        print(n_acc_1.unit)
        print(n_p.unit)
        raise Exception("***WARNING*** n_acc and n_p must have same units")

    idx = find_nearest(ene,ene_2[0])

    n_fin = s**(2/3) * np.concatenate((n_acc_1.value, np.maximum(n_p.value[idx:],n_acc_2.value)))
    
    return n_fin
