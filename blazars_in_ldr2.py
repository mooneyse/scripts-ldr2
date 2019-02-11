#!/usr/bin/env python3

'''Get the blazars from BZCAT that matched with the LDR2 data that has been
analysed to date.'''

import argparse
import glob
import pandas as pd

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '08 February 2019'

def blazars_in_ldr2(files):
    '''From the filenames of the images Tim sent on, read the names of the
    sources.'''

    blazars = glob.glob(files)  # get list of filenames
    print('Searching directory {}.'.format(files))
    all = []

    for blazar in blazars:
        blazar = blazar.split('/')[-1]  # first split
        blazar = blazar.split('-P')[0]  # second split
        all.append(blazar)  # get list of blazar names

    return all


def read_in_bzcat(catalogue, blazars, output):
    '''Read in the BZCAT.'''

    pd.set_option('expand_frame_repr', False)
    pd.set_option('display.max_columns', None)

    df = pd.read_csv(catalogue, sep='|')  # read in bzcat

    for i in range(len(blazars)):
        blazars[i] = ' ' + blazars[i] + ' '  # need to buffer by a space

    db = df[df[' Source name '].isin(blazars)]  # filter bzcat by ldr2 blazars
    db.to_csv(output)

    return db


def sample_properties(data):
    '''Print some of the properties of the sample.'''

    source_types = [' BL Lac ', ' QSO RLoud flat radio sp. ']
    valid_redshift = data[~data[' Redshift '].str.contains('\?')]  # remove uncertain redshift values
    valid_redshift = valid_redshift[valid_redshift[' Redshift '].astype('float64') > 0]  # remove zeroes
    gamma_detected = data[data[' Fermi flux1-100 GeV(ph/cm2/s)'] > 0]

    statistics = {'bl_lac': data[data[' Source classification '] == ' BL Lac '].shape[0], # number of bl lacs
                  'fsrq': data[data[' Source classification '] == ' QSO RLoud flat radio sp. '].shape[0], # number of fsrqs
                  'other': data[~data[' Source classification '].isin(source_types)].shape[0],
                  'flux_mean': data[' Flux density 1.4/0.843GHz(mJy) '].mean(),  # average flux density
                  'flux_minimum': data[' Flux density 1.4/0.843GHz(mJy) '].min(),  # average flux density
                  'redshift': valid_redshift[' Redshift '].astype('float64').mean(),  # average redshift
                  'gamma': gamma_detected.shape[0],  # gamma-ray detected sources
                  'total': data.shape[0]}  # total source count

    print('Of the {total} sources, there are {bl_lac} BL Lacs, {fsrq} FSRQs, and {other} other sources.'.format(**statistics))
    print('The average 1.4 GHz flux density is {flux_mean:.0f} mJy.'.format(**statistics))
    print('The minimum 1.4 GHz flux density is {flux_minimum:.0f} mJy.'.format(**statistics))
    print('The sensitivity of LDR2 is 0.07 mJy and we expect the indices to be flat.')
    print('The average redshift is {redshift:.2f}.'.format(**statistics))
    print('There are {gamma} gamma-ray 3FGL sources.'.format(**statistics))


def main():
    '''Get the LDR2 blazar names from the images Tim produced.'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-f',
                        '--files',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/images/png/*.png',
                        help='Files containing the data')

    parser.add_argument('-c',
                        '--catalogue',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/catalogues/bzcat.txt',
                        help='BZCAT as a text file')

    parser.add_argument('-o',
                        '--output',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/catalogues/ldr2-bzcat.csv',
                        help='File to save the BZCAT filtered by LDR2.')

    args = parser.parse_args()
    files = args.files
    catalogue = args.catalogue
    output = args.output

    blazars = blazars_in_ldr2(files)
    data = read_in_bzcat(catalogue, blazars, output)
    sample_properties(data)

if __name__ == '__main__':
    main()
