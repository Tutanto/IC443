import numpy as np
import astropy.units as u
from astropy.table import QTable
from astropy.io import ascii
from reacc import reacc

def make_table(s,co,br):

    # Prepare an energy array for saving the particle distribution
    energy = np.logspace(-6, 5, 100) * u.GeV

    n_pr_reacc = reacc(energy,co,br,s) / u.cm**3/u.eV

    n_pr_reacc = np.delete(n_pr_reacc,0)
    energy = np.delete(energy,0)

    name1 = 'energy'
    name2 = 'flux'

    data = QTable([energy, n_pr_reacc], names=(name1, name2))
    ascii.write(data, 'proton_reacc.dat', format='ipac')

    return

