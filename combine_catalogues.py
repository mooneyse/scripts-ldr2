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

def clean(x):
    '''Get the BZCAT source name from the FITS filename.'''
    return x.split('-P')[0]


def combine_catalogues(fluxes, redshifts, polarisation, fermi):
    '''Combine data into a new CSV using Pandas.'''

    # TODO Crossmatch 4FGL
    
    pd.set_option('expand_frame_repr', False)
    pd.set_option('display.max_columns', None)

    df_polarisation = pd.read_csv(polarisation)
    df_polarisation.rename(columns={'Source_name': 'BZCAT name'}, inplace=True)  # rename column for merging
    df_flux = pd.DataFrame(fluxes, columns=['BZCAT name', 'Total flux 144 MHz (mJy)', 'Peak flux 144 MHz (mJy)'])
    df_redshift = pd.DataFrame(redshifts, columns=['NED name', 'BZCAT name', 'Redshift'])

    df_flux['BZCAT name'] = df_flux['BZCAT name'].apply(clean)  # remove path from name
    df_flux_redshift = pd.merge(df_flux, df_redshift, on='BZCAT name', how='inner')
    df = pd.merge(df_flux_redshift, df_polarisation, on='BZCAT name', how='left')

    drop = ['Source_Name', '# col1', 'id', 'Redshift_y', 'Separation_2',
            'ID_ra', 'ID_dec', 'Separation', 'RA',  'DEC', 'RAdeg', 'DECdeg',
            'CatID', 'TargetID', 'z_best_source', 'z_best', 'Peak_flux',
            'Total_flux', 'LOFAR?']

    for column in drop:  # remove duplicate columns
        df.drop(column, axis=1, inplace=True)  # axis=0 for rows, axis=1 for columns

    df.rename(columns={'RA_(J2000.0)': 'RA',
                       'Dec_(J2000.0)': 'Dec',
                       'Redshift_x': 'Redshift',
                       'Flux_density_1.4/0.843GHz(mJy)': 'Flux density 1.4 GHz (mJy)',
                       'Flux_density143GHz(mJy)': 'Flux density 14 3GHz (mJy)',
                       'X-ray_flux0.1-2.4_keV(1.e-12_cgs)': 'Flux 1.3 keV (1e-12cgs)',
                       'Fermi_flux1-100_GeV(ph/cm2/s)': 'Flux 51 GeV (ph/cm2/s)',
                       'Source_classification': 'Classification'},
                       inplace=True)

    print(df)


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
