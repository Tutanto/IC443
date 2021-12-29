import numpy as np
import astropy.units as u
from astropy.table import QTable
from astropy.io import ascii
from reacc_kinetic import reacc_kinetic

def make_table(s,co_p,br_p,co_e,br_e,path):

    # Prepare an energy array for saving the particle distribution
    energy = np.logspace(-6, 5, 100) * u.GeV

    n_pr_reacc = reacc_kinetic(energy,co_p,br_p,s,particle = 'proton') / u.cm**3/u.eV
    n_el_reacc = reacc_kinetic(energy,co_e,br_e,s,particle = 'electron') / u.cm**3/u.eV

    n_pr_reacc = np.delete(n_pr_reacc,0)
    n_el_reacc = np.delete(n_el_reacc,0)
    energy = np.delete(energy,0)

    name1 = 'energy'
    name2 = 'flux'

    data = QTable([energy, n_pr_reacc], names=(name1, name2))
    ascii.write(data, path+'/proton_reacc.dat', format='ipac', overwrite='True')

    data = QTable([energy, n_el_reacc], names=(name1, name2))
    ascii.write(data, path+'/electron_reacc.dat', format='ipac', overwrite='True')

    return

