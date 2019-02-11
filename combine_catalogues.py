#!/usr/bin/env python3

'''Combine the total flux and peak flux I derived from the FITS files, the
redshifts I got from NED, and the polarisation data with the BZCAT sources
using Pandas.'''

import argparse
import pandas as pd
from flux_from_fits import flux_from_fits
from query_ned import query_ned

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '11 February 2019'

# TODO add 4FGL crossmatches

def combine_catalogues(fluxes, redshifts, polarisation, fermi):
    '''Combine data into a new CSV using Pandas.'''

    df_pol = pd.read_csv(polarisation)
    df_flux = pd.DataFrame(fluxes, columns=['BZCAT name', 'Total flux', 'Peak flux'])
    df_redshift = pd.DataFrame(redshifts, columns=['NED name', 'BZCAT name', 'Redshift'])

    # name = flux_name.split('-P')[0]

    print(df_flux.head)
    print(df_redshift.head)


def main():
    '''Combine the total flux and peak flux I derived from the FITS files, the
    redshifts I got from NED, and the polarisation data with the BZCAT sources
    using Pandas.'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-f',
                        '--fits',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/images/fits/*.fits',
                        help='FITS files containing the data')

    parser.add_argument('-c',
                        '--catalogue',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/catalogues/ldr2-bzcat.csv',
                        help='CSV of BZCAT sources in LDR2')

    parser.add_argument('-p',
                        '--polarisation',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/catalogues/ldr2-bzcat.csv',
                        help='CSV of polarised LDR2 sources in BZCAT')

    parser.add_argument('-F',
                        '--fermi',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/catalogues/gll_psc_8year_v6.fit',
                        help='4FGL')

    args = parser.parse_args()
    fits = args.fits
    catalogue = args.catalogue
    polarisation = args.polarisation
    fermi = args.fermi

    fluxes = flux_from_fits(fits)
    redshifts = query_ned(catalogue)
    combine_catalogues(fluxes, redshifts, polarisation, fermi)


if __name__ == '__main__':
    main()
