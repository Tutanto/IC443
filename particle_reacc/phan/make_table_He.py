import numpy as np
import astropy.units as u
from astropy.table import QTable
from astropy.io import ascii
from reacc_He import reacc_He

def make_table_He(s,co_p,br_p):

    # Prepare an energy array for saving the particle distribution
    energy = np.logspace(-6, 5, 100) * u.GeV

    n_he_reacc = reacc_He(energy,co_p,br_p,s) / u.cm**3/u.eV

    n_he_reacc = np.delete(n_he_reacc,0)
    energy = np.delete(energy,0)

    name1 = 'energy'
    name2 = 'flux'

    data = QTable([energy, n_he_reacc], names=(name1, name2))
    ascii.write(data, 'he_reacc.dat', format='ipac')

    return

