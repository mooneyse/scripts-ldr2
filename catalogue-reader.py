#!/usr/bin/env python3

'''Print information from a catalogue that is in FITS format.'''

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '07 February 2019'

def catalogue_reader(fits):
    table = Table.read(fits, format='fits')
    df = table.to_pandas()
    df.to_csv(fits[:-4] + 'csv')

    pd.set_option('expand_frame_repr', False)
    pd.set_option('display.max_columns', None)
    print(df)

    # ds = df[df['S_Code']=='S']
    # print('Sources fit with a single Guassian:', ds.shape[0])
    # ds = df[df['S_Code']!='S']
    # print('Sources fit with more than a single Guassian:', ds.shape[0])


def main():
    '''Print information from a catalogue that is in FITS format.'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-f',
                        '--fits',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/catalogues/VanEck_HETDEX_LoTSS_76_final.fits',
                        help='FITS file containing the data')

    args = parser.parse_args()
    fits = args.fits

    catalogue_reader(fits)


if __name__ == '__main__':
    main()
