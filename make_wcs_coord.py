'''
Read a 2D .fits map and create HEADER information to add WCS coordinates.
A .reg file is needed, with a region of the same shape and size of the map.
Creates a new map with the new HEADER.
'''
import pyregion
import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.io import fits

image = 'maps/sum_flux_map_reacc2_mcl+acl_ott.fits'  # path to the map
region = 'regions/big_box.reg'  # path to region

hdu = fits.open(image)[0]  # open the map
r = pyregion.open(region).as_imagecoord(hdu.header)  # open the region files

dl = 1.6890676e+17 * u.cm  # physical dimension of pixel
dist = 1.5*u.kpc  # distance from IC443
# angular diameter of the pixel
delta = (2*np.arctan(dl/(2*dist.to(u.cm)))).to(u.deg)

# (207,352,330) => explosion centre simulation domain coordinates
center = [207, 330]

# center and size of region
x_cntr = r[0].coord_list[0]
y_cntr = r[0].coord_list[1]
x_size = r[0].coord_list[2]
y_size = r[0].coord_list[3]

index_x1 = int(x_cntr - x_size/2)
index_x2 = int(x_cntr + x_size/2)
index_y1 = int(y_cntr - y_size/2)
index_y2 = int(y_cntr + y_size/2)

# creates header
wcs_input_dict = {
    'CTYPE1': 'RA---TAN',
    'CUNIT1': 'deg',
    'CDELT1': -delta.value,
    'CRPIX1': center[0]-index_x1,
    'CRVAL1': 94.32,
    'CTYPE2': 'DEC--TAN',
    'CUNIT2': 'deg',
    'CDELT2': delta.value,
    'CRPIX2': center[1]-index_y1,
    'CRVAL2': 22.38,
}
wcs_helix_dict = WCS(wcs_input_dict)

for key, value in wcs_input_dict.items():
    hdu.header.set(key, value)

# write a new map with the header information
hdu.writeto('maps/sum_flux_map_reacc2_mcl+acl_ott_wcs.fits', clobber=True)
