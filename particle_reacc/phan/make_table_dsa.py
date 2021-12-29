import numpy as np
import astropy.units as u
from astropy.table import QTable
from astropy.io import ascii
from dsa import dsa

def make_table_dsa(A,norm,a,co,K):
    # Prepare an energy array for saving the particle distribution
    energy = np.logspace(-4, 4, 100) * u.GeV

    pr = dsa(energy,A,norm,a,co,particle = 'proton') /u.eV
    el = dsa(energy,A*K,norm,a,co,particle = 'electron') /u.eV

    name1 = 'energy'
    name2 = 'flux'

    data = QTable([energy, pr], names=(name1, name2))
    ascii.write(data, 'proton_dsa.dat', format='ipac')

    data = QTable([energy, el], names=(name1, name2))
    ascii.write(data, 'electron_dsa.dat', format='ipac')

    return

