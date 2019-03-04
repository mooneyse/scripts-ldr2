#!/usr/bin/env python2.7

'''Use PyBDSF (http://www.astron.nl/citt/pybdsf/) to make the noise maps from a
list of fits files by fitting 2D Gaussians and removing them from the data.'''

import argparse
import bdsf
import glob
import numpy as np
import matplotlib.pyplot as plt
# import pandas as pd
# from astropy import units as u
# from astropy.coordinates import SkyCoord

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '04 March 2019'

def get_noise(fits, directory='/mnt/closet/ldr2-blazars/images'):
    '''Get the noise of an image.'''

    name = fits.split('/')[-1].split('P')[0][:-1]
    skymodel = directory + '/skymodel/' + name + '.skymodel'
    noise_map = directory + '/noise/' + name + '.fits'
    image = bdsf.process_image(fits)
    image.write_catalog(outfile=skymodel, bbs_patches='source', format='bbs', srcroot=name)
    image.export_image(outfile=noise_map, img_type='gaus_resid')
    return name, np.mean(image.rms_arr)


def make_noise_maps(directory):
    '''Use PyBDSF to fit the sources.'''

    fits_files = glob.glob(directory + '/*')
    rms = []
    for fits in fits_files:
        rms.append(get_noise(fits))
        break
    print(rms)
    # TODO write noise to file and match it with the csv data
    # what I want to do is say if the max(abs) > ~0.6 then flag=BAD, otherwise flag=GOOD
    # and save this with the noise value.


def main():
    '''Use PyBDSF (http://www.astron.nl/citt/pybdsf/) to make the noise maps
    from a list of fits files by fitting 2D Gaussians and removing them from
    the data.'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-d', '--directory', type=str, help='Directory of FITS files', default='/mnt/closet/ldr2-blazars/images/fits')

    args = parser.parse_args()
    directory = args.directory

    make_noise_maps(directory)


if __name__ == '__main__':
    main()
