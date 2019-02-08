#!/usr/bin/env python3

'''Query NED for the spectral information for a source, given the position.'''

import argparse
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt
from astroquery.ned import Ned
from astropy import coordinates

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '08 February 2019'

def plot_sed(frequency, flux):
    '''Plot the SED from NED data.'''

    plt.loglog(frequency, flux, marker = '.', ls='None', color='black')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Flux density (Jy)')
    plt.show()


def query_ned(catalogue):
    '''Query NED given the source position.'''

    df = pd.read_csv(catalogue, sep=',')  # read in ldr2 bzcat sources

    for ra, dec in zip(df[' RA (J2000.0) '], df[' Dec (J2000.0) ']):
        coordinate = coordinates.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
        result_table = Ned.query_region(coordinate, radius=10 * u.arcsec, equinox='J2000.0')
        object_name = str(result_table[0]['Object Name'])[2:-1]

        photometry = Ned.get_table(object_name, table='photometry')
        frequency = photometry['Frequency']  # hertz
        flux = photometry['Flux Density']  # jansky
        plot_sed(frequency, flux)
        break

    return result_table


def main():
    '''Query NED for the spectral information for a source, given the position.
    '''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-c',
                        '--catalogue',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/catalogues/ldr2-bzcat.csv',
                        help='Files containing the data')

    args = parser.parse_args()
    catalogue = args.catalogue

    query_ned(catalogue)


if __name__ == '__main__':
    main()
