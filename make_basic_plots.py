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


def g(x, pos):
    f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    return r'${}$'.format(f._formatSciNotation('%1.10e' % x))


def my_plot(my_directory='/mnt/closet/ldr2'):
    """Make some plots.

    Parameters
    ----------
    my_directory : string
        Working directory.
    """

    df = pd.read_csv(f'{my_directory}/catalogues/final.csv')
    # f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    # g = lambda x, pos: '${}$'.format(f._formatSciNotation('%1.10e' % x))

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
    plt.errorbar(df['Extent (kpc)'], df['LDR2-to-NVSS index'],
                 xerr=df['Extent (kpc)'] * 0.2,
                 yerr=df['LDR2-to-NVSS index'] * 0.2,
                 marker='o', ls='none', color='k')
    plt.xlabel('Extent (kpc)', fontsize=20)
    plt.ylabel(r'$\alpha_{\mathrm{MHz-GHz}}$', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    # plt.show()
    plt.savefig(f'{my_directory}/images/extent-v-index.png')
    plt.clf()

    # plot extent (kpc) against radio luminosity
    plt.errorbar(df['Extent (kpc)'], df['Luminosity (W/Hz)'], marker='o',
                 xerr=df['Extent (kpc)'] * 0.2,
                 yerr=df['Luminosity (W/Hz)'] * 0.2,
                 ls='none', color='k')
    plt.xlabel('Extent (kpc)', fontsize=20)
    plt.ylabel(r'$L$ (W Hz$^{-1}$)', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks([0, 1e26, 2e26, 3e26, 4e26, 5e26],
               ['0', r'$1 \times 10^{26}$', r'$2 \times 10^{26}$',
                r'$3 \times 10^{26}$', r'$4 \times 10^{26}$',
                r'$5 \times 10^{26}$'], fontsize=20)
    plt.tight_layout()
    # plt.show()
    plt.savefig(f'{my_directory}/images/extent-v-luminosity.png')
    plt.clf()

    # plot luminosity against redshift
    plt.errorbar(df['Redshift'], df['Luminosity (W/Hz)'], marker='o',
                 xerr=df['Redshift'] * 0.2,
                 yerr=df['Luminosity (W/Hz)'] * 0.2,
                 ls='none', color='k')
    plt.xlabel('Redshift', fontsize=20)
    plt.ylabel(r'$L$ (W Hz$^{-1}$)', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks([0, 1e26, 2e26, 3e26, 4e26, 5e26],
               ['0', r'$1 \times 10^{26}$', r'$2 \times 10^{26}$',
                r'$3 \times 10^{26}$', r'$4 \times 10^{26}$',
                r'$5 \times 10^{26}$'], fontsize=20)
    plt.tight_layout()
    # plt.show()
    plt.savefig(f'{my_directory}/images/redshift-v-luminosity.png')
    plt.clf()

    # plot luminosity against spectral index
    plt.errorbar(df['LDR2-to-NVSS index'], df['Luminosity (W/Hz)'], marker='o',
                 xerr=df['LDR2-to-NVSS index'] * 0.2,
                 yerr=df['Luminosity (W/Hz)'] * 0.2,
                 ls='none', color='k')
    plt.xlabel(r'$\alpha_{\mathrm{MHz-GHz}}$', fontsize=20)
    plt.ylabel(r'$L$ (W Hz$^{-1}$)', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.yticks([0, 1e26, 2e26, 3e26, 4e26, 5e26],
               ['0', r'$1 \times 10^{26}$', r'$2 \times 10^{26}$',
                r'$3 \times 10^{26}$', r'$4 \times 10^{26}$',
                r'$5 \times 10^{26}$'], fontsize=20)
    plt.tight_layout()
    # plt.show()
    plt.savefig(f'{my_directory}/images/index-v-luminosity.png')


def main():
    """Make some plots.
    """
    my_plot()


if __name__ == '__main__':
    main()
