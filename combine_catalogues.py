#!/usr/bin/env python3

'''Combine the total flux and peak flux I derived from the FITS files, the
redshifts I got from NED, and the polarisation data with the BZCAT sources
using Pandas.'''

import argparse
from flux_from_fits import flux_from_fits
from query_ned import query_ned


__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '11 February 2019'

# TODO add 4FGL crossmatches




def main():
    '''Combine the total flux and peak flux I derived from the FITS files, the
    redshifts I got from NED, and the polarisation data with the BZCAT sources
    using Pandas.'''

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

    fluxes = flux_from_fits(files)
    # for flux in fluxes:
    #     print(flux.split('-P')[0])
    for peak, total, name in fluxes:
        print(name.split('-P')[0])


if __name__ == '__main__':
    main()
