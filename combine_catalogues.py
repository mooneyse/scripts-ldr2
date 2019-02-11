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

def clean(x):
    '''Get the BZCAT source name from the FITS filename.'''

    return x.split('-P')[0]


def combine_catalogues(fluxes, redshifts, polarisation, fermi):
    '''Combine data into a new CSV using Pandas.'''

    pd.set_option('expand_frame_repr', False)
    pd.set_option('display.max_columns', None)

    df_polarisation = pd.read_csv(polarisation)
    df_polarisation.rename(columns={'Source_name': 'BZCAT name'}, inplace=True)  # rename column for merging
    df_flux = pd.DataFrame(fluxes, columns=['BZCAT name', 'Flux density (144 MHz)', 'Peak flux (144 MHz)'])
    df_redshift = pd.DataFrame(redshifts, columns=['NED name', 'BZCAT name', 'Redshift'])

    df_flux['BZCAT name'] = df_flux['BZCAT name'].apply(clean)  # remove path from name
    df_flux_redshift = pd.merge(df_flux, df_redshift, on='BZCAT name', how='inner')
    df = pd.merge(df_flux_redshift, df_polarisation, on='BZCAT name', how='left')

    for column in ['Source_Name', '# col1', 'id', 'Redshift_y', 'Separation_2', 'Separation', 'RA',  'DEC']:  # remove duplicate columns
        df.drop(column, axis=1, inplace=True)  # axis=0 for rows, axis=1 for columns

    print(list(df))
    print(df.head(5))


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
                        default='/mnt/closet/ldr2-blazars/catalogues/ldr2-bzcat-vaneck.csv',
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
