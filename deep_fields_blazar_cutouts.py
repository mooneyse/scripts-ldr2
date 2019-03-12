#!/usr/bin/env python3

'''Get images of the blazars in the Deep Fields data.'''

import matplotlib as mpl
mpl.use('Agg')

import argparse
import aplpy
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy import coordinates
from astropy import units as u

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '12 March 2019'

def make_cut_out_image(sources, field, radius=1 / 60, cmap='viridis', vmin=0,
                       output='/data5/sean/deep-fields'):
    '''Make a cut-out image of a given source.'''

    df = pd.read_csv(sources)
    # print(df.head)

    for name, ra, dec, peak_flux in zip(df['name'], df['RA'], df['DEC'], df['Peak_flux']):
        image = aplpy.FITSFigure(field)
        image.show_colorscale(cmap=cmap, vmin=vmin, vmax=peak_flux)
        image.recenter(ra, dec, radius=radius)
        image.add_colorbar()
        image.set_title(name)
        image.save(output + '/' + name + '.png')
        sys.exit()


def main():
    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=formatter_class)

    parser.add_argument('-s', '--sources', required=False, type=str,
                        help='CSV of the BZCAT',
                        default='/data5/sean/deep-fields/bootes-lockman-hole-blazars.csv')

    parser.add_argument('-b', '--bootes', required=False, type=str,
                        help='FITS image of the Bootes field',
                        default='/data5/sean/deep-fields/bootes/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked.scaled.fits')

    parser.add_argument('-l', '--lockman', required=False, type=str,
                        help='FITS image of the Lockman Hole',
                        default='/data5/sean/deep-fields/lockman-hole/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked.scaled.fits')

    args = parser.parse_args()
    sources = args.sources
    bootes = args.bootes
    lockman = args.lockman

    fields = [bootes, lockman]  # elais n1

    for field in fields:
            make_cut_out_image(sources, field)


if __name__ == '__main__':
    main()
