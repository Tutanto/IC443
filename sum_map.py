
from astropy.io import fits

image = 'maps/flux_map_reacc2_mcl_ott.fits'  # path to first map
hdu = fits.open(image)[0]  # open the map file
data1 = hdu.data  # store the data

image = 'maps/flux_map_reacc2_acl_ott.fits'
hdu = fits.open(image)[0]
data2 = hdu.data

# overwrite the data of the second map with the sum of the 2 maps
hdu.data = data1 + data2
# save the sum in a new file using the same header of the second map
hdu.writeto('maps/sum_flux_map_reacc2_mcl+acl_ott.fits', clobber=True)
