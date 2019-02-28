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
    # good (P), marginal (S), not enough correlated flux (X), bad station (-)
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

    '''Results:
    ID       RA            Dec           Date        Time      Flux           Code  RA (deg)    Dec (deg)  Separation
    L618324  08:54:48.870  20:06:30.701  2017-10-05  06:59:24  PPPPPPPXPPXX-  56    133.703625  20.108528  0.000000   OJ 287
    L618322  08:53:02.750  20:04:21.601  2017-10-05  06:59:24  XXXXXXXXXXXXX  56    133.261458  20.072667  0.416807
    L618058  08:53:13.400  19:30:52.999  2017-10-05  07:27:32  PPPPPPPXPPXX-  36    133.305833  19.514722  0.701901   Yes
    L618340  08:56:39.800  20:57:44.302  2017-10-05  06:59:24  XXXXXXXXXXXX-  56    134.165833  20.962306  0.957225
    L618344  08:56:57.210  21:11:44.300  2017-10-05  06:59:24  PPPPPPPXPPXX-  56    134.238375  21.195639  1.196740   Yes
    L618064  08:57:31.570  18:58:38.201  2017-10-05  07:27:32  XXXXXXXXXXXX-  36    134.381542  18.977278  1.299174
    L618238  08:49:25.950  18:40:00.199  2017-10-05  07:06:26  SSSSXSXXPXXX-  50    132.358125  18.666722  1.920821
    L618326  08:46:17.030  20:46:29.701  2017-10-05  06:59:24  XXXXXXXXXXXX-  56    131.570958  20.774917  2.106521
    L618346  08:50:37.800  22:06:15.700  2017-10-05  06:59:24  XXXXXXXXXXXX-  56    132.657500  22.104361  2.221637
    L618300  09:05:13.170  20:29:50.302  2017-10-05  06:59:24  PPSPPPPXPPXPX  56    136.304875  20.497306  2.470390
    '''


if __name__ == '__main__':
    main()
