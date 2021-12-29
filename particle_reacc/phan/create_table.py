import numpy as np
import astropy.units as u
from astropy.table import QTable
from astropy.io import ascii
from reacc_kinetic import reacc_kinetic

# Prepare an energy array for saving the particle distribution
energy = np.logspace(-3, 8, 200) * u.MeV
s = 12.4
co = 14 * u.GeV
br = 85 * u.GeV

n_pr_reacc = reacc_kinetic(energy,co.to(u.MeV),br.to(u.MeV),s,particle = 'proton') / u.cm**3/u.eV
n_el_reacc = reacc_kinetic(energy,co.to(u.MeV),br.to(u.MeV),s,particle = 'electron') / u.cm**3/u.eV

name1 = 'energy'
name2 = 'flux'

data = QTable([energy, n_pr_reacc], names=(name1, name2))
ascii.write(data, 'proton_reacc.dat', format='ipac')

data = QTable([energy, n_el_reacc], names=(name1, name2))
ascii.write(data, 'electron_reacc.dat', format='ipac')

