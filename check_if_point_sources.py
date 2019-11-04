#!/usr/bin/env python3

"""Make some plots.
"""

import pandas as pd

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '06 November 2019'


def is_compact(S_int, S_peak, RMS):
    if (S_int / S_peak) < 1.25 + 3.1 * (S_peak / RMS) ** -0.53:
        return True
    else:
        return False


def check_if_point_source(my_directory='/mnt/closet/ldr2'):
    # read in the total flux, the peak flux, and the noise

    # check the pulsars first
    psr = pd.read_csv(f'{my_directory}/catalogues/pulsars-10asec-match.csv')
    for name, S_int, S_peak, RMS in zip(psr['NAME'], psr['Total_flux'],
                                        psr['Peak_flux'], psr['Isl_rms']):
        print(f'{name}: {is_compact(S_int=S_int, S_peak=S_peak, RMS=RMS)}')

    print()

    blc = pd.read_csv(f'{my_directory}/catalogues/final-20-with-point-sources.csv')
    for name, S_int, S_peak, RMS, S_code in zip(blc['Source name'],
                                                blc['Total_flux'],
                                                blc['Peak_flux'],
                                                blc['Isl_rms'],
                                                blc['S_Code']):
        print(f'{name}\t{is_compact(S_int=S_int, S_peak=S_peak, RMS=RMS)}\t'
              f'{S_code}')


def main():
    """Check if sources are point sources according to the equation given in
    the caption of Fig. 7 of Shimwell et al. (2019).
    """
    check_if_point_source()


if __name__ == '__main__':
    main()
