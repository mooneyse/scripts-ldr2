#!/usr/bin/env python3

'''Plot the LDR2 radio contours over the SDSS optical images which are fetched
from the website with an API.'''

import warnings
warnings.filterwarnings('ignore')  # supress warnings

import argparse
import aplpy
import glob
import pandas as pd
import matplotlib.pyplot as plt
from astroquery.sdss import SDSS
from astropy import coordinates as coords

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '16 February 2019'

def query_sdss(ra, dec, unit='deg', frame='icrs', band='i', spectro=False):
    '''Query SDSS. Many filters are available. See this URL for more:
    http://skyserver.sdss.org/dr2/en/proj/advanced/color/sdssfilters.asp.'''

    position = coords.SkyCoord(ra, dec, unit=unit, frame=frame)
    id = SDSS.query_region(position, spectro=spectro)
    images = SDSS.get_images(matches=id, band=band)

    return images[0]  # use first image


def make_plot(fits_file, df, format='png', north=True, figsize=12, radius=120,
              cmap='bone_r', vmin=0, vmax=1, radio_cmap='YlOrRd_r',
              radio_levels=[4, 8, 16, 32],
              save_directory='/mnt/closet/ldr2-blazars/images/contour'):
    '''Plot contours from one FITS file over the image of another. Radius is
    given in arcseconds.'''

    # house-keeping
    blazar_name = fits_file.split('/')[-1].split('-P')[0]
    blazar_row = df[df['name'] == blazar_name]
    ra, dec = float(blazar_row['ra']), float(blazar_row['dec'])
    save = '{}/{}.{}'.format(save_directory, blazar_name, format)
    try:
        image_file = query_sdss(ra, dec)
    except:
        print('SDSS match for {} not found.'.format(blazar_name))
        return

    # make image
    image = aplpy.FITSFigure(image_file, north=north, figsize=(figsize, figsize))
    image.show_colorscale(cmap=cmap, vmin=vmin, vmax=vmax)
    image.recenter(ra, dec, radius=radius / (60 * 60))  # degrees
    image.add_scalebar(radius / (60 * 60 * 2))
    image.scalebar.set_label(str(int(radius / 2)) + '"')
    image.show_grid()
    image.set_title(blazar_name)
    image.show_contour(fits_file, cmap=radio_cmap, levels=radio_levels)  # mJy
    image.save(save)
    print('Image saved at {}.'.format(save))


def radio_over_optical(radio_directory, bzcat_csv):
    '''Get the list of sources and for each, call a function to plot the radio
    contours over the optical image.'''

    df = pd.read_csv(bzcat_csv)  # get catalogue which has positions we need
    df.rename(columns={' Source name ': 'name', ' RA (J2000.0) ': 'ra',
                       ' Dec (J2000.0) ': 'dec'}, inplace=True)
    df['name'] = df['name'].str.strip()  # remove spaces from blazar names
    fits_files = glob.glob('{}/*'.format(radio_directory))  # get fits files

    for fits_file in fits_files:
        make_plot(fits_file, df)


def main():
    '''Plot the LDR2 radio contours over the SDSS optical images which are
    fetched from the website with an API.'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-d', '--directory', type=str, help='BZCAT CSV',
                        default='/mnt/closet/ldr2-blazars/images/fits')

    parser.add_argument('-c', '--csv', type=str, help='BZCAT CSV',
                        default='/mnt/closet/ldr2-blazars/catalogues/ldr2-bzcat.csv')

    args = parser.parse_args()
    directory = args.directory
    csv = args.csv

    radio_over_optical(radio_directory=directory, bzcat_csv=csv)


if __name__ == '__main__':
    main()
