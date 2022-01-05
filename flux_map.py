'''
This script reads initial and final density data cubes from PLUTO and computes
the gamma-ray flux in all computational cells inside a region defined by a .reg file.
Then integrate along y-axes.
The script can work on multiple processors
'''
import os
import sys
import pyregion
import numpy as np
import astropy.units as u
from scipy import integrate
import astropy.io.fits as fits
from astropy.table import QTable
from multiprocessing import Pool
from naima.models import PionDecay, InverseCompton, Bremsstrahlung, TableModel


def enumerate2D(array1, array2):
    '''
    Parameters
    ----------
    array1 : numpy.array
        Generic array.
    array2 : numpy.array
        Generic array, same dimension of array1.

    Yields
    ------
    indexes :
        list of indexes of the array
    data :
        Content of arrays

    This function is the same as numpy.ndenumerate, but works
    on two different arrays in the same time

    '''
    assert array1.shape == array2.shape, "Error - dimensions."
    for indexes, data in np.ndenumerate(array1):
        yield indexes, data, array2[indexes]


def make_gamma(frange):
    '''
    Parameters
    ----------
    frange : list
        first and last x-index for the i-th processor

    Returns
    -------
    int_flux : a flux map


    This function loops over the initial and final density cubes and
    computes che flux map

    '''
    path = 'tables/particles_' + \
        str(frange[0])+'-'+str(frange[1]
                               )  # directory to store protons&electrons tables
    os.makedirs(path, exist_ok=True)  # succeeds even if directory exists.
    # empty 3D array to fill
    var = np.zeros((z_f-z_i, int(y_size), int(x_size/n_processes)))

    for indexes, data1, data2 in enumerate2D(dens_ini[z_i:z_f, index_y1:index_y2, frange[0]:frange[1]], dens_f[z_i:z_f, index_y1:index_y2, frange[0]:frange[1]]):

        if (data1 > 0 and data2 > 0):  # loop only in non-empy cells

            s = (data2/data1)/4  # compute compression factor

            # create particle distribution table
            make_table(s, cut_p, p_break, cut_e, e_break, path)

            # Read protons and electrons energy distribution tables
            spectrum_prot = QTable.read(
                path+'/proton_reacc.dat', format='ipac')
            spectrum_elec = QTable.read(
                path+'/electron_reacc.dat', format='ipac')

            hadr = TableModel(
                energy=spectrum_prot['energy'],
                values=spectrum_prot['flux'],
                amplitude=vol,
            )

            lept = TableModel(
                energy=spectrum_elec['energy'],
                values=spectrum_elec['flux'],
                amplitude=vol,
            )
            ###

            nh = data2 * u.cm**-3  # cell density

            # Set emission mechanisms
            PP = PionDecay(hadr, nh=nh)
            IC = InverseCompton(lept, seed_photon_fields=["CMB"])
            BRE = Bremsstrahlung(lept, n0=nh)

            # Compute flux emitted from cell in energy range = energy
            pion = PP.flux(energy, distance=distance)
            inverse = IC.flux(energy, distance=distance)
            bre = BRE.flux(energy, distance=distance)
            flux = pion + inverse + bre  # total flux from all emission mechanisms

            # integrate flux in the energy range using Simpson method
            int_flux = integrate.simps(flux.value, energy.value)
            var[indexes] = int_flux  # store integral result in var

    # collapse 3D array in 2D along "line of sight"
    flux_map_section = var.sum(1)

    return flux_map_section


if __name__ == '__main__':

    # Global variables

    distance = 1.5 * u.kpc  # distance of IC443
    energy = np.logspace(6, 14, 100) * u.eV  # flux map energy range

    vol = 4.826809e+51 * u.cm**3  # single computational cell volume
    cut_p = 5000 * u.GeV  # cut off energy protons distribution
    cut_e = 5000 * u.GeV  # cut off energy electrons distribution
    p_break = 18 * u.GeV  # proton break energy
    e_break = 40 * u.GeV  # electron break energy

    n_processes = 50  # number of cpus to perform the calculation
    z_i = 139  # starting z_index
    z_f = 599  # ending z_index
    #####

    # path to set of functions to compute re-accelerated particles
    path = "./particle_reacc/phan"
    sys.path.insert(1, path)  # adding path to sys
    # import function to create re-accelerated particles energy distribution table
    from make_table import make_table

    # Load data cubes
    image = 'original_data/IC443_V5/density_sh4_mcl_0f.fits'  # path to initial data cube
    hdu = fits.open(image)[0]
    dens_ini = hdu.data  # initial density data cube
    header = hdu.header  # header info
    size_cube = len(dens_ini)  # data cube pixel size

    image = 'original_data/IC443_V5/density_sh4_mclf.fits'  # path to final data cube
    hdu = fits.open(image)[0]
    dens_f = hdu.data  # final density data cube

    # Load .reg to select a portion of the simulated domain where compute the gamma-ray map
    r = pyregion.open('regions/big_box.reg').as_imagecoord(header)

    # Extract center and size (in pixel) of the .reg region
    x_cntr = r[0].coord_list[0]
    y_cntr = r[0].coord_list[1]
    x_size = r[0].coord_list[2]
    y_size = r[0].coord_list[3]

    # Indexes where the .reg region starts & ends in the simulated domain
    index_x1 = int(x_cntr - x_size/2)
    index_x2 = int(x_cntr + x_size/2)
    index_y1 = int(y_cntr - y_size/2)
    index_y2 = int(y_cntr + y_size/2)

    # Set the "parallel" proces.
    # Divide the whole simulated domain inside the .reg region in n_processes stripes of x-width = sub_len.
    # Each stripe is handled by a single process
    pool = Pool(processes=n_processes)
    sub_len = int(x_size / n_processes)
    d = {}
    # Create an object containg starting & ending x-index for each stripe
    for x in range(1, n_processes+1):
        if x == n_processes:
            d["DATA_SET_{0}".format(x)] = [index_x1 + sub_len*(x-1), index_x2]
        else:
            d["DATA_SET_{0}".format(x)] = [index_x1 +
                                           sub_len*(x-1), index_x1 + sub_len*x]

    # Run make_gamma function in n_processes cpus each one with different frange value (defined in d.values)
    res = pool.map(make_gamma, list(d.values()))

    # concatenate output from each cpu
    final_flux = np.concatenate(res, axis=1)
    hdu.data = final_flux
    # save result in .fits with same header of dens_f
    hdu.writeto('test.fits', clobber=True)
