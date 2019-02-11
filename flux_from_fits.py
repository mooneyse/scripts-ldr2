#!/usr/bin/env python3

'''Sum the total flux at the centre of a FITS file.'''

import argparse
import glob
import numpy as np
from astropy.io import fits

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '11 February 2019'

def box(x, c):
    '''Dimensions of the centre box.'''
    return int(x / 2 - c / 2), int(x / 2 + c / 2)


def flux_from_fits(files, cutout=100):
    '''Sum the total flux at the centre of a FITS file.'''

    blazars = glob.glob(files)  # get list of filenames
    fluxes = []
    i = 0
    for blazar_fits in blazars:
        hdu_list = fits.open(blazar_fits)
        data = np.squeeze(hdu_list[0].data)  # remove single-dimensional entries
        l, L = box(data.shape[0], cutout)  # length
        w, W = box(data.shape[1], cutout)  # width
        printing = {'file': blazar_fits.split('/')[-1],
                    'peak': np.max(data[l:L, w:W]),
                    'total': np.sum(data[l:L, w:W])}
        hdu_list.close()
        print('File: {file}, Peak flux: {peak} mJy/beam, Total flux: {total} mJy/beam'.format(**printing))
        fluxes.append(list(printing.values()))
        i += 1
        if i > 3:
            break
    return fluxes


def main():
    '''Sum the total flux at the centre of a FITS file.'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-f',
                        '--files',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/images/fits/*.fits',
                        help='Files containing the data')

    args = parser.parse_args()
    files = args.files

    flux_from_fits(files)


if __name__ == '__main__':
    main()
