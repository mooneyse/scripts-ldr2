#!/usr/bin/env python2.7

'''Use PyBDSF (http://www.astron.nl/citt/pybdsf/) to make the noise maps from a
list of fits files by fitting 2D Gaussians and removing them from the data.'''

import argparse
import bdsf
import csv
import glob
import numpy as np
import matplotlib.pyplot as plt
# import pandas as pd
# from astropy import units as u
# from astropy.coordinates import SkyCoord

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '04 March 2019'

def get_noise(fits, directory):
    '''Get the noise of an image.'''

    name = fits.split('/')[-1].split('P')[0][:-1]
    skymodel = directory + '/images/skymodel/' + name + '.skymodel'
    noise_map = directory + '/images/noise/' + name + '.fits'

    image_max, image_rms, residual_max = float('nan'), float('nan'), float('nan')

    try:
        image = bdsf.process_image(fits)
        # print image.__dict__.keys()  # print object attributes
        image.write_catalog(outfile=skymodel, bbs_patches='source', format='bbs', srcroot=name)
        image.export_image(outfile=noise_map, img_type='gaus_resid')
        # using max instead of np.max below gives an error
        image_max = np.max(np.abs(image.image_arr))
        image_rms = np.mean(image.rms_arr)
        residual_max = np.max(np.abs(image.resid_gaus_arr))
    except:
        print 'Failed for %s.' % name

    return name, image_max, image_rms, residual_max


def make_noise_maps(directory):
    '''Use PyBDSF to fit the sources.'''

    fits_files = glob.glob(directory + '/images/fits/*.fits')
    rms = []
    for fits in fits_files:
        rms.append(get_noise(fits, directory))

    csv_file = directory + '/catalogues/noise.csv'
    with open(csv_file, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['Name', 'Max_image', 'RMS', 'Max_residual'])
        writer.writerows(rms)


def main():
    '''Use PyBDSF (http://www.astron.nl/citt/pybdsf/) to make the noise maps
    from a list of fits files by fitting 2D Gaussians and removing them from
    the data.'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-d', '--directory', type=str, help='Working directory', default='/mnt/closet/ldr2-blazars')

    args = parser.parse_args()
    directory = args.directory

    make_noise_maps(directory)


if __name__ == '__main__':
    main()
