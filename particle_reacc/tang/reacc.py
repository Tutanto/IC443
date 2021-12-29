import numpy as np
import astropy.units as u
from astropy.constants import c, m_p, m_e
from scipy import integrate
from ambient import ambient_CR

def reacc(ene,co,br,s):

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
    
    E_p = m_p * c**2
    rest = E_p.to('GeV')
    rest = rest.value
        
   # Parameters protons
    const = 1.56e-10/u.cm**3/u.GeV
    a = 1.03
    b = 1.21
    d = -3.18

    def beta(kin):
        factor = (kin+rest)/rest
        beta = factor/np.sqrt(1+(factor)**2)
        return beta

    def k2p(kin):
        return (np.sqrt(kin**2 + 2*rest*kin))#/cc

    f = lambda x: (k2p(x))**(2*(alpha-1)) * (x + rest) * (k2p(x)**a)/beta(k2p(x))**3 * ((k2p(x)**a + 0.77**b)/(1 + 0.77**b))**d

    np.vectorize(f)

    def fint(a,z):
        return integrate.quad(f, a, z)

    fint=np.vectorize(fint)

    an_1=fint(min(ene.value), ene_1)
    an_2=fint(min(ene.value), ene_2)

    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    n_acc_1 = an_1[0] * (alpha+2) * (k2p(ene_1))**(-alpha) * np.exp(-(k2p(ene_1))/co) * const
    n_acc_2 = an_2[0] * (alpha+2) * (k2p(ene_2))**(-alpha) * br/(k2p(ene_2)) * np.exp(-(k2p(ene_2))/co) * const 
    
    n_acc_1 = n_acc_1.to(u.eV**(-1)/u.cm**3)
    n_acc_2 = n_acc_2.to(u.eV**(-1)/u.cm**3)

    n_p = ambient_CR(ene)
    n_p = n_p.to(u.eV**(-1)/u.cm**3)

    # Check if n_acc as same units as n_p
    if not n_acc_1.unit == n_p.unit or not n_acc_2.unit == n_p.unit:
        print(n_acc_1.unit)
        print(n_p.unit)
        raise Exception("***WARNING*** n_acc and n_p must have same units")

    idx = find_nearest(ene,ene_2[0])

    n_fin = s**(2/3) * np.concatenate((n_acc_1.value, np.maximum(n_p.value[idx:],n_acc_2.value)))
    return n_fin
