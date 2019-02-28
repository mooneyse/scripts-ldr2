#!/usr/bin/env python3

'''Search LBCS catalogue for sources nearest a given position.'''

import argparse
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '28 February 2019'

def lbcs_search(ra, dec, lbcs, nearest=5, sep='  ', engine='python'):
    '''Search LBCS.'''

    pd.set_option('expand_frame_repr', False)
    pd.set_option('display.max_columns', None)

    names = ['ID', 'RA', 'Dec', 'Date', 'Time', 'Flux', 'Code', 'RA (deg)', 'Dec (deg)']
    df = pd.read_csv(lbcs, names=names, sep=sep, engine=engine)  # get catalogue which has positions we need
    separations = []
    c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')

    for ra_, dec_ in zip(df['RA (deg)'], df['Dec (deg)']):
        c_ = SkyCoord(ra=ra_ * u.degree, dec=dec_ * u.degree, frame='icrs')
        separations.append(c.separation(c_).degree)

    df['Separation'] = separations
    df.sort_values('Separation', inplace=True)
    print(df.head(10))


def main():
    '''Search LBCS catalogue for sources nearest a given position.'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-r', '--ra', type=float, help='Right ascension to search',
                        default=133.703625)

    parser.add_argument('-d', '--dec', type=float, help='Declination to search',
                        default=20.108528)

    parser.add_argument('-l', '--lbcs', type=str, help='LBCS',
                        default='/home/sean/Downloads/workbooks/lbcs_stats.sum')

    args = parser.parse_args()
    ra = args.ra
    dec = args.dec
    lbcs = args.lbcs

    lbcs_search(ra=ra, dec=dec, lbcs=lbcs)


if __name__ == '__main__':
    main()
