#!/usr/bin/env python3

'''Query NED for the spectral information for a source, given the position.'''

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy import coordinates
from astropy import units as u
from PIL import Image, ImageDraw
from astroquery.ned import Ned as ned

from catalogue_reader import catalogue_reader

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '11 February 2019'

def flux_from_fits(catalogue, savefig='/mnt/closet/ldr2-blazars/images/sed'):
    '''Query NED given the source position.'''

    pass


def main():
    '''Query NED for the spectral information for a source, given the position.
    '''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-f',
                        '--fits',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/catalogues/ldr2-bzcat.csv',
                        help='Files containing the data')

    args = parser.parse_args()
    fits = args.fits

    flux_from_fits(fits)


if __name__ == '__main__':
    main()
