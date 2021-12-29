import numpy as np
import astropy.units as u
from astropy.constants import c, m_p, m_e
from scipy.integrate import quad

def dsa(ene,A,norm,a,co,particle):

    # Check if input energy is GeV
    if not ene.unit is u.GeV or not co.unit is u.GeV or not norm.unit is u.GeV:
        raise Exception("Please use GeV unit")
    # Check if
#    if not A.unit is 1/u.eV:
#        raise Exception("Please use 1/eV unit")

    ene = ene.value
    norm = norm.value
    co = co.value
    A = A.value
    cc = c.to(u.cm/u.s)
#    cc = cc.value

    if particle == 'proton':
        E_p = m_p * c**2
        rest = E_p.to('GeV')
        rest = rest.value

    if particle == 'electron':
        E_e = m_e * c**2
        rest = E_e.to('GeV')
        rest = rest.value

    co_k = np.sqrt(rest**2 + co**2) - rest

    def k2p(ene):
        p = (np.sqrt(ene**2 + 2*rest*ene))#/cc.value
        return p
    def dp_dk(ene):
#        var = (ene + rest)/(cc.value*np.sqrt(ene*(ene+2*rest)))
        var = (ene + rest)/(np.sqrt(ene*(ene+2*rest)))
        return var
    def cutpwl(x,norm,A,a,co_k):
        dn_dk = A * (k2p(x)/norm)**(-a) * np.exp(-(k2p(x))/co_k) * dp_dk(x)
#       dn_dk_2 = A * (k2p(ene_2)/norm) **(-b) * (br/norm)**(b-a) * np.exp(-(k2p(ene_2))/co) * dp_dk(ene_2)
        return dn_dk

    I = quad(cutpwl, 1, 1e7, args=(norm,A,a,co_k))

    print(I)

    fl = cutpwl(ene,norm,A,a,co_k)

    return fl

