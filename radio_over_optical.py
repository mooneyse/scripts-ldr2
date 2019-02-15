#!/usr/bin/env python3

import warnings
warnings.filterwarnings('ignore')  # supress warnings

import aplpy
from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astropy.io import fits
from astropy.wcs import WCS
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# TODO add this all into a function

radio_directory = '/mnt/closet/ldr2-blazars/images/fits'
save_directory = '/mnt/closet/ldr2-blazars/images/contour'
bzcat_csv = '/mnt/closet/ldr2-blazars/catalogues/ldr2-bzcat.csv'
band = 'i'  # i = near infrared
cmap = 'bone_r'
vmin, vmax = 0, 1
format = 'png'
radio_cmap = 'YlOrRd_r'
radio_levels = [4, 8, 16, 32]  # mJy
radius = 120  # arcseconds
figsize = 12
unit = 'deg'
frame = 'icrs'

df = pd.read_csv(bzcat_csv)  # get catalogue which has positions we need
df.rename(columns={' Source name ': 'name', ' RA (J2000.0) ': 'ra',
                   ' Dec (J2000.0) ': 'dec'}, inplace=True)
df['name'] = df['name'].str.strip()  # remove spaces from blazar names
fits_files = glob.glob('{}/*'.format(radio_directory))  # get the fits files

for fits_file in fits_files:
    blazar_name = fits_file.split('/')[-1].split('-P')[0]

    blazar_row = df[df['name'] == blazar_name]
    ra, dec = float(blazar_row['ra']), float(blazar_row['dec'])

    save = '{}/{}.{}'.format(save_directory, blazar_name, format)

    # http://skyserver.sdss.org/dr2/en/proj/advanced/color/sdssfilters.asp
    position = coords.SkyCoord(ra, dec, unit=unit, frame=frame)
    id = SDSS.query_region(position, spectro=False)
    image_file = SDSS.get_images(matches=id, band=band)[0]  # use first image

    image = aplpy.FITSFigure(image_file, north=True, figsize=(figsize, figsize))
    image.show_colorscale(cmap=cmap, vmin=vmin, vmax=vmax)
    image.recenter(ra, dec, radius=radius / (60 * 60))  # degrees
    image.add_scalebar(radius / (60 * 60 * 2))
    image.scalebar.set_label(str(int(radius / 2)) + '"')
    image.show_grid()
    image.set_title(blazar_name)
    print(fits_file)
    image.show_contour(fits_file, cmap=radio_cmap, levels=radio_levels)
    image.save(save)
    break
