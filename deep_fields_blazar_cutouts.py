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

def make_cut_out_image(sources, radius=1 / 60, cmap='viridis', vmin=0,
                       output='/data5/sean/deep-fields/images',
                       contours=[5, 10, 20]):
    '''Make a cut-out image of a given source.'''

    df = pd.read_csv(sources)

    for source_name, ra, dec, peak_flux, field in zip(df['Source name'],
        df['RA'], df['DEC'], df['Peak_flux'], df['Field']):
        # TODO add noise field to include contours

        if field == 'Bootes':
            field_file = '/data5/sean/deep-fields/bootes/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked.scaled.fits'
        elif field == 'Lockman Hole':
            field_file = '/data5/sean/deep-fields/lockman-hole/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked.scaled.fits'

        image = aplpy.FITSFigure(field_file)
        image.recenter(ra, dec, radius=radius)
        image.show_colorscale(cmap=cmap, vmin=vmin, vmax=peak_flux,
                              stretch='arcsinh')
        # image.show_contour(levels=[contour * noise for contour in contours],
        #                    colors='white', overlap=True)
        image.add_colorbar()
        # image.add_scalebar(radius / 4)
        # image.scalebar.set_label('15"')
        image.set_title(source_name)
        image.save(output + '/' + source_name + '.png')


def main():
    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-s', '--sources', required=False, type=str,
                        help='CSV of the BZCAT',
                        default='/data5/sean/deep-fields/bootes-lockman-hole-blazars.csv')

    args = parser.parse_args()
    sources = args.sources

    make_cut_out_image(sources)


if __name__ == '__main__':
    main()
