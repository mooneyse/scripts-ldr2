#!/usr/bin/env python3

'''Plot postage stamp images of LDR2 BL Lacs.'''

import warnings
warnings.filterwarnings('ignore')  # supress warnings

import aplpy
import argparse
import math
import os
import sys
import numpy as np
import pandas as pd
import astropy.io
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw

__author__ = 'Sean Mooney'
__date__ = '27 June 2019'

def files_in_directory(directory):
    '''Get the files in a given directory.'''
    files = [directory + '/' + f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
    return files

def get_levels(fits, sigma=[3, 5]):
    '''Get the the noise for the FITS image by taking the standard deviation of the residual.'''
    residual_fits = f'/mnt/closet/ldr2-blazars/images/noise/{fits[41:55]}.fits'
    # using the residuals is not ideal so the noise is over estimated
    data = np.squeeze(astropy.io.fits.open(residual_fits)[0].data)
    return np.std(data), [s * np.std(data) for s in sigma]  # mJy


def make_plot(fits, has_z, z='', z_flag='', new_diameter='', vmin=0):
    '''Get a png image from a FITS file. The size is the width of the square in arcseconds.'''
    data = np.squeeze(astropy.io.fits.open(fits)[0].data)
    pixels, _ = data.shape  # they are square
    diameter = 12  # arcminutes

    if not new_diameter:
        new_diameter = diameter

    # find cutout region
    scaling = new_diameter / diameter
    new_pixels = pixels * scaling
    start = np.rint((pixels / 2) - (new_pixels / 2)).astype(int)
    end = np.rint((pixels / 2) + (new_pixels / 2)).astype(int)
    cutout = data[start:end, start:end]

    # get contour levels
    noise, levels = get_levels(fits)

    # make figure
    plt.figure(figsize=(15, 12))
    plt.imshow(cutout, origin='lower', vmin=0, vmax=np.max(cutout))
    cbar = plt.colorbar()
    plt.contour(cutout, levels=levels, colors=['white', 'magenta'])
    plt.title(f'{os.path.basename(fits)[:14]}, z = {z_flag} {z}', fontsize=20)
    plt.xlabel(f'Radius = {new_diameter / 2}"; RMS = {noise:.3f} mJy', fontsize=20)
    plt.tick_params(axis='both', which='both', left=False, right=False, labelleft=False, bottom=False, top=False, labelbottom=False)
    cbar.set_label(f'Jy beam\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}')
    plt.tight_layout()
    extra_dir = 'with-z' if has_z else 'without-z'
    ending = 'z' if has_z else 'no-z'
    save = f'/mnt/closet/ldr2-blazars/images/png/bzb-for-fra/{extra_dir}/{os.path.basename(fits)[:14]}-{ending}.png'
    plt.savefig(save)
    plt.clf()
    print(f'Imaged saved at {save}.')


def get_z(blazar, z_csv, no_z_csv):

    blazar = blazar[:4] + ' ' + blazar[4:]  # have to add a space
    df_z = pd.read_csv(z_csv)
    df_no_z = pd.read_csv(no_z_csv)

    if blazar in list(df_z['Name']):
        row = df_z[df_z['Name'] == blazar]
        z = row['Redshift'].values[0]
        flag = ''
        has_z = True
        print(f'{blazar} has a good z (z = {z}).')

    elif blazar in list(df_no_z['Name']):
        row = df_no_z[df_no_z['Name'] == blazar]
        z = row['Redshift'].values[0]
        flag = row['Redshift flag'].values[0]
        has_z = False
        if math.isnan(z):
            z = ''
            print(f'{blazar} does not have a z.')
        if math.isnan(flag):
            flag = ''
        else:
            print(f'{blazar} does not have a good z (z = {flag} {z}).')

    else:
        z = ''
        flag = ''
        has_z = False
        print('BL Lac not in either list.')

    return has_z, z, flag


def main():
    '''Plot postage stamp images of LDR2 BL Lacs.'''
    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=formatter_class)
    parser.add_argument('-d', '--directory', required=False, type=str, default='/mnt/closet/ldr2-blazars/images/fits/bzb', help='FITS files to be plotted.')
    parser.add_argument('-z', '--z_csv', required=False, type=str, default='/mnt/closet/ldr2-blazars/catalogues/bl_lacs_with_z.csv', help='CSV of BL Lacs with redshift.')
    parser.add_argument('-n', '--no_z_csv', required=False, type=str, default='/mnt/closet/ldr2-blazars/catalogues/bl_lacs_without_z.csv', help='CSV of BL Lacs without redshift.')
    args = parser.parse_args()
    directory = args.directory
    z_csv = args.z_csv
    no_z_csv = args.no_z_csv

    fits_list = files_in_directory(directory=directory)

    for fits in fits_list:
        try:
            has_z, z, z_flag = get_z(os.path.basename(fits)[:14], z_csv, no_z_csv)
            make_plot(fits=fits, has_z=has_z, new_diameter=2, z=z, z_flag=z_flag)
        except:
            print(f'Failed for {fits}.')


if __name__ == '__main__':
    main()
