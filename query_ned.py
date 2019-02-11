#!/usr/bin/env python3

'''Query NED for the spectral information for a source, given the position.'''

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy import coordinates
from astropy import units as u
from PIL import Image, ImageDraw
from astroquery.ned import Ned as ned

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '08 February 2019'

def TeV_to_MHz(TeV=50):
    '''Use astropy to find the equivalency between MHz and TeV, which I used to
    find the ideal frequency axis range for the blazar SEDs.'''

    exact = (TeV * u.TeV).to(u.MHz, equivalencies=u.spectral()).value
    return 10 ** np.round(np.log10(exact))  # nearest power of 10


def plot_sed(source, frequency, flux, savefig, fontsize=16):
    '''Plot the SED from NED data.'''

    nu_f_nu = np.array(frequency) * np.array(flux)

    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True
    mpl.rcParams['xtick.major.size'] = 10
    mpl.rcParams['xtick.major.width'] = 2
    mpl.rcParams['xtick.minor.size'] = 5
    mpl.rcParams['xtick.minor.width'] = 2
    mpl.rcParams['ytick.major.size'] = 10
    mpl.rcParams['ytick.major.width'] = 2
    mpl.rcParams['ytick.minor.size'] = 5
    mpl.rcParams['ytick.minor.width'] = 2
    mpl.rcParams['axes.linewidth'] = 2

    plt.figure(figsize=(12, 8))
    plt.loglog(frequency, nu_f_nu, marker='.', ls='None', color='black')
    plt.xlim(1e7, TeV_to_MHz())  # 10 Hz to 50 TeV
    plt.xlabel('Frequency (Hz)', fontsize=fontsize)#, labelpad=100)
    plt.ylabel(r'$\nu \cdot f_{\nu}$ (Jy Hz)', fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.minorticks_off()
    plt.title(source, fontsize=fontsize)

    plt.tight_layout()
    fname = '{}/{}.png'.format(savefig, source)  # leading and trailing spaces
    plt.savefig(fname, bbox_inches='tight')
    plt.close()  # close open figures
    print('SED saved at {}.'.format(fname))


def query_ned(catalogue, savefig='/mnt/closet/ldr2-blazars/images/sed'):
    '''Query NED given the source position.'''

    df = pd.read_csv(catalogue, sep=',')  # read in ldr2 bzcat sources

    for source, ra, dec in zip(df[' Source name '], df[' RA (J2000.0) '], df[' Dec (J2000.0) ']):
        source = source.strip()  # remove leading and trailing spaces
        coordinate = coordinates.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
        result_table = ned.query_region(coordinate, radius=10 * u.arcsec, equinox='J2000.0')
        object_name = str(result_table[0]['Object Name'])[2:-1]
        redshift = result_table[0]['Redshift']

        try:  # query ned for photometry information
            photometry = ned.get_table(object_name, table='photometry')
            frequency = photometry['Frequency']  # hertz
            flux = photometry['Flux Density']  # jansky
            plot_sed(source, frequency, flux, savefig=savefig)
        except:  # if there is no photometry, create a blank image
            string = 'No photometry on NED for {}'.format(source)
            print(string + '.')
            length, width = 1189, 790  # the size of the matplotlib figures
            blank = Image.new('RGB', (length, width), color='white')
            draw = ImageDraw.Draw(blank)
            draw.text((length / 2, width / 2), string, fill='black')  # approximate centre
            blank.save('{}/{}.png'.format(savefig, source))


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
