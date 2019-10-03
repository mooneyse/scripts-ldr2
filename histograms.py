#!/usr/bin/env python3

"""Make some plots.
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '06 September 2019'


def my_plot(my_directory='/mnt/closet/ldr2'):
    """Make some plots.

    Parameters
    ----------
    my_directory : string
        Working directory.
    """

    df = pd.read_csv(f'{my_directory}/catalogues/final.csv')
    f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    g = lambda x, pos: '${}$'.format(f._formatSciNotation('%1.10e' % x))

    plt.figure(figsize=(13.92, 8.60)).patch.set_facecolor('white')
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    mpl.rcParams['xtick.major.size'] = 10
    mpl.rcParams['xtick.major.width'] = 2
    mpl.rcParams['xtick.minor.size'] = 5
    mpl.rcParams['xtick.minor.width'] = 2
    mpl.rcParams['ytick.major.size'] = 10
    mpl.rcParams['ytick.major.width'] = 2
    mpl.rcParams['ytick.minor.size'] = 5
    mpl.rcParams['ytick.minor.width'] = 2
    mpl.rcParams['axes.linewidth'] = 2

    # plot extent (kpc) against spectral index
    plt.plot(df['LDR2-to-NVSS index'], df['Extent (kpc)'], marker='o',
             ls='none', color='k')
    plt.xlabel(r'$\alpha_{\mathrm{MHz-GHz}}$', fontsize=20)
    plt.ylabel('Extent (kpc)', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    # plt.show()
    plt.savefig(f'{my_directory}/images/fig1-alpha-v-extent.png')
    plt.clf()

    # plot extent (kpc) against radio luminosity
    plt.plot(df['Luminosity (W/Hz)'], df['Extent (kpc)'], marker='o',
             ls='none', color='k')
    plt.xlabel(r'$L$ (W Hz$^{-1}$)', fontsize=20)
    plt.ylabel('Extent (kpc)', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.gca().xaxis.set_major_formatter(mticker.FuncFormatter(g))
    # plt.show()
    plt.savefig(f'{my_directory}/images/fig2-luminosity-v-extent.png')
    plt.clf()

    # plot luminosity against redshift
    plt.plot(df['Luminosity (W/Hz)'], df['Redshift'], marker='o',
             ls='none', color='k')
    plt.xlabel(r'$L$ (W Hz$^{-1}$)', fontsize=20)
    plt.ylabel('Redshift', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.gca().xaxis.set_major_formatter(mticker.FuncFormatter(g))
    # plt.show()
    plt.savefig(f'{my_directory}/images/fig3-luminosity-v-redshift.png')
    plt.clf()

    # plot luminosity against spectral index
    plt.plot(df['Luminosity (W/Hz)'], df['LDR2-to-NVSS index'], marker='o',
             ls='none', color='k')
    plt.xlabel(r'$L$ (W Hz$^{-1}$)', fontsize=20)
    plt.ylabel(r'$\alpha_{\mathrm{MHz-GHz}}$', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.gca().xaxis.set_major_formatter(mticker.FuncFormatter(g))
    # plt.show()
    plt.savefig(f'{my_directory}/images/fig4-luminosity-v-alpha.png')


def main():
    """Make some plots.
    """
    my_plot()


if __name__ == '__main__':
    main()
