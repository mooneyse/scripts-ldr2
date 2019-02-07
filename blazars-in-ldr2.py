#!/usr/bin/env python3

'''Download the catalogues that contain LDR2 blazars that Tim imaged for me on
01 February 2019 when in Dublin. It queries the lofar-surveys.org website and
gets the catalogue.'''

import argparse
import glob
import subprocess

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '07 February 2019'

def get_ldr2_catalogues(files):
    '''From the filenames of the images Tim sent on, read the pointings so I can
    download the catalogues from the website.'''

    blazars = glob.glob(files)  # get list of filenames
    print('Searching directory {}.'.format(files))
    all_mosaics, urls = [], []

    for blazar in blazars:
        all_mosaics.append(blazar.split('-')[2])  # get list of mosaics

    mosaics = sorted(list(set(all_mosaics)))  # get unique ordered list
    print('{} blazars found in {} mosaics.'.format(len(blazars), len(mosaics)))

    for mosaic in mosaics:  # build the urls from the mosaics
        url = 'https://lofar-surveys.org/downloads/DR2/mosaics/{}/mosaic.cat.fits'.format(mosaic)
        urls.append([mosaic, url])

    for m, u in urls:
        command = ['wget', '--user', 'surveys', '--password', '150megahertz', '--output-document', m + '.fits', u]

        try:
            output = subprocess.check_output(command)
            print('{}.fits downloaded successfully from {}.'.format(m, u))
        except:
            print('{}.fits from {} failed.'.format(m, u))


def main():
    '''Download the catalogues that contain LDR2 blazars.'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-f',
                        '--files',
                        required=False,
                        type=str,
                        default='/mnt/closet/blazars-shimwell/5BZ*.png',
                        help='FITS file containing the data')

    args = parser.parse_args()
    files = args.files

    get_ldr2_catalogues(files)


if __name__ == '__main__':
    main()
