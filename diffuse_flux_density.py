#!/usr/bin/env python3

'''Measure the flux density for a source at different levels above the
noise.'''

import aplpy
import argparse
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs import WCS, utils
from astropy import coordinates
from astropy import units as u

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '17 March 2019'

def coordinates_for_blazar(blazar_name, csv):
    '''Get the position of a source from a CSV file.'''

    df = pd.read_csv(csv)
    df.set_index('Source name', inplace=True)
    ra = df.loc[blazar_name, 'RA']
    dec = df.loc[blazar_name, 'DEC']

    return ra, dec

def shift_and_crop(ra, dec, fits_file, size=60):
    '''Shift the FITS file so it centres on the source and crop it around this
    point.'''

    hdu = fits.open(fits_file)[0]  # load the image
    position = SkyCoord(ra, dec, unit='deg', frame='icrs')
    wcs = WCS(hdu.header)  # load the wcs
    pixel = utils.skycoord_to_pixel(position, wcs=wcs)
    data = np.squeeze(hdu.data)
    cutout = Cutout2D(data, position=pixel, size=size)#, wcs=wcs)
    cutout_data = cutout.data
    hdu.data = cutout_data  # put the cutout image in the fits hdu
    # hdu.header.update(cutout.wcs.to_header())  # update header with cutout WCS
    fits_split = fits_file.split('.')
    cutout_file = fits_split[0] + '-cutout.' + fits_split[-1]
    hdu.writeto(cutout_file, clobber=True)
    print('Cut-out written to {}.'.format(cutout_file))

    return cutout_data


def plot_it(cutout_data, rms_noise, diffuse_sigma, core_sigma):
    '''Plot the FITS cut-out file.'''

    vmax = np.max(cutout_data)
    diffuse_level = [diffuse_sigma * rms_noise]
    core_level = [core_sigma * rms_noise]

    plt.figure()
    plt.imshow(cutout_data, cmap='viridis', vmin=0, vmax=vmax, origin='lower')
    plt.contour(cutout_data, levels=diffuse_level, colors='white', alpha=0.5)
    plt.contour(cutout_data, levels=core_level, colors='red', alpha=0.5)
    plt.show()


def measure_flux_density(cutout_data, minimum_sigma, rms_noise):
    '''Measure the flux density.'''

    return cutout_data[cutout_data > minimum_sigma * rms_noise].sum()


def give_results(diffuse_and_core, core_only):
    '''Format the results.'''

    diffuse_only = diffuse_and_core - core_only
    print('The core and jet flux density is {:.3f} Jy/beam.'.format(core_only))
    print('The diffuse flux density is {:.3f} Jy/beam.'.format(diffuse_only))


def main():
    '''Measure the flux density for a source at different levels above the
    noise.'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-f',
                        '--fits_file',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/deep-fields/bootes-image.fits',
                        help='FITS file containing the data')

    parser.add_argument('-c',
                        '--csv',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/deep-fields/bootes-lockman-hole-blazars.csv',
                        help='CSV file with the blazar data')

    parser.add_argument('-b',
                        '--blazar_name',
                        required=False,
                        type=str,
                        default='5BZQJ1437+3519',
                        help='Name of the blazar')

    parser.add_argument('-D',
                        '--diffuse_sigma',
                        required=False,
                        type=float,
                        default=5,
                        help='Minimum sigma value for the diffuse emission')

    parser.add_argument('-C',
                        '--core_sigma',
                        required=False,
                        type=float,
                        default=320,
                        help='Minimum sigma value for the core emission')

    parser.add_argument('-r',
                        '--rms_noise',
                        required=False,
                        type=float,
                        default=0.000136,
                        help='RMS noise in the FITS file')

    args = parser.parse_args()
    fits_file = args.fits_file
    csv = args.csv
    blazar_name = args.blazar_name
    diffuse_sigma = args.diffuse_sigma
    core_sigma = args.core_sigma
    rms_noise = args.rms_noise

    ra, dec = coordinates_for_blazar(blazar_name=blazar_name, csv=csv)
    cutout_data = cropped_fits = shift_and_crop(ra=ra, dec=dec,
                                                fits_file=fits_file)
    plot_it(cutout_data=cutout_data, rms_noise=rms_noise,
            diffuse_sigma=diffuse_sigma, core_sigma=core_sigma)
    diffuse_and_core = measure_flux_density(cutout_data=cutout_data,
                                            minimum_sigma=diffuse_sigma,
                                            rms_noise=rms_noise)  # jansky
    core_only = measure_flux_density(cutout_data=cutout_data,
                                            minimum_sigma=core_sigma,
                                            rms_noise=rms_noise)  # jansky
    give_results(diffuse_and_core=diffuse_and_core, core_only=core_only)


if __name__ == '__main__':
    main()
