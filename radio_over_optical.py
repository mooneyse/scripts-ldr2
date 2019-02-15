#!/usr/bin/env python3

import warnings
warnings.filterwarnings('ignore')  # supress warnings

import aplpy
from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

pos = coords.SkyCoord(187.277915, 2.052388, unit='deg', frame='icrs')
xid = SDSS.query_region(pos)#, spectro=False)
# print(xid)

# http://skyserver.sdss.org/dr2/en/proj/advanced/color/sdssfilters.asp
im = SDSS.get_images(matches=xid, band='i')  # i = near infrared
# print(im[0][0])
# # image_data = fits.getdata(im[0][0], ext=0)
# # print(im[0].info())#.info)
# data = im[0][0][1].data
# image_data = im[0][0]
# print(type(image_data))
#
# print(image_data)
hdul = im[0]
# print(hdul[0].header['ALT'])
hdr = hdul[0].header

data = hdul[0].data
            # print(data.shape)
            # print(data.dtype.name)
            # plt.imshow(data, cmap='gray')
            # plt.show()
print(repr(hdr))
print(np.min(data), np.max(data))

print(hdr['CRPIX1'], hdr['CRPIX2'])  # x, y of ref pixel
print(hdr['CRVAL1'], hdr['CRVAL2'])  # ra, dec of ref pixel

# image = aplpy.FITSFigure(data, convention='calabretta')  # , figsize=(16, 16)
# image = aplpy.FITSFigure('/home/sean/Downloads/fits/sdss-test.fits', hdu=0)#, figsize=(16, 16))
# w1,w2= image.pixel2world(100,100)
# print(image.data)



# print(w1,w2)
# image.set_xaxis_coord_type('scalar')
# image.recenter(187.277915, 2.052388, radius=60 / 60 / 60)
# image.show_colorscale(cmap='viridis', vmin=0, vmax=1)
# image.add_grid()
# plt.show()






# image = aplpy.FITSFigure(fits)  # north=True
# image.recenter(centre[0], centre[1], radius=radius / 60 / 60)
# image.show_colorscale(cmap=cmap, vmin=vmin, vmax=vmax)
# image.frame.set_linewidth(0)
# image.hide_axis_labels()
# image.hide_tick_labels()
#
# save = os.path.splitext(fits)[0] + '.' + format
# image.save(save, dpi=dpi, max_dpi=max_dpi)


# im[0].close()

# plt.imshow(image_data, cmap='gray')

image_file = '/home/sean/Downloads/fits/sdss-test.fits'
# hdu_list = fits.open(image_file)
# hdu_list.info()
# image_data = hdu_list[0].data
#
#
# print(type(image_data))
# print(image_data.shape)
#
# plt.imshow(image_data, cmap='gray')
# plt.colorbar()
# plt.show()

from astropy.wcs import WCS
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)

# plt.subplot(projection=wcs)
# plt.imshow(hdu.data,  vmin=-0, vmax=1, origin='lower')
# plt.grid(color='white', ls='solid')
# plt.xlabel('Galactic Longitude')
# plt.ylabel('Galactic Latitude')
# plt.show()


ax = plt.subplot(projection=wcs)

ax.imshow(hdu.data, vmin=-0, vmax=1, origin='lower')

ax.coords.grid(True, color='white', ls='solid')
# ax.coords[0].set_axislabel('Galactic Longitude')
# ax.coords[1].set_axislabel('Galactic Latitude')

# overlay = ax.get_coords_overlay('fk5')
# overlay.grid(color='white', ls='dotted')
# overlay[0].set_axislabel('Right Ascension (J2000)')
# overlay[1].set_axislabel('Declination (J2000)')
plt.show()
