'''
This script collapses a 3D data cube in a 2D map along "line-of-sight"
'''

import astropy.io.fits as fits  

image = 'original_data/density_sh4_mcl_0.fits'
hdu = fits.open(image)[0]  #ccube
hdu.data = hdu.data.sum(1)
hdu.writeto('maps/dens_map_sh4_mcl_0.fits', clobber=True)
