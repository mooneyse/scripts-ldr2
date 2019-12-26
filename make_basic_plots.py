#!/usr/bin/env python3

"""Make some plots.
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from math import sqrt, sin, exp
from matplotlib.lines import Line2D

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '06 September 2019'


def get_dl_and_kpc_per_asec(z, H0=70, WM=0.26, WV=0.74):
    """Get luminosity distance. See
    """
    WR = 0.  # Omega(radiation)
    WK = 0.  # Omega curvaturve = 1-Omega(total)
    c = 299792.458  # velocity of light in km/sec
    DTT = 0.5  # time from z to now in units of 1/H0
    age = 0.5  # age of Universe in units of 1/H0
    zage = 0.1  # age of Universe at redshift z in units of 1/H0
    DCMR = 0.0  # comoving radial distance in units of c/H0
    DA = 0.0  # angular size distance
    DA_Mpc = 0.0
    kpc_DA = 0.0
    DL = 0.0  # luminosity distance
    DL_Mpc = 0.0
    a = 1.0  # 1/(1+z), the scale factor of the Universe
    az = 0.5  # 1/(1+z(object))
    h = H0 / 100.
    WR = 4.165E-5 / (h * h)  # includes 3 massless neutrino species,
    WK = 1 - WM - WR - WV
    az = 1.0 / (1 + 1.0 * z)
    age = 0.
    n = 1000  # number of points in integrals
    for i in range(n):
        a = az * (i + 0.5) / n
        adot = sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
        age = age + 1. / adot
    zage = az * age / n
    DTT = 0.0
    DCMR = 0.0
    for i in range(n):
        a = az + (1 - az) * (i + 0.5) / n
        adot = sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
        DTT = DTT + 1. / adot
        DCMR = DCMR + 1. / (a * adot)
    DTT = (1. - az) * DTT / n
    DCMR = (1. - az) * DCMR / n
    age = DTT + zage
    ratio = 1.00
    x = sqrt(abs(WK)) * DCMR
    if x > 0.1:
        if WK > 0:
            ratio = 0.5 * (exp(x) - exp(-x)) / x
        else:
            ratio = sin(x) / x
    else:
        y = x * x
        if WK < 0:
            y = -y
        ratio = 1. + y / 6. + y * y / 120.
    return (c / H0) * ((az * ratio * DCMR) / (az * az)) * 3.086e22, ((c / H0) * az * ratio * DCMR) / 206.264806

import numpy as np
redshift,lum=[],[]
for i in range(0, 1000):
    z = i / 1000
    dist,s = get_dl_and_kpc_per_asec(z=i/1000)
    if z == 0.64:
        print(dist,s)
    L = (0.07 *5 * 0.001 * 1E-26 * 4. * np.pi * dist * dist)/((1. + z) ** (1))
    redshift.append(z)
    lum.append(L)
lum07=[]
for i in range(0, 1000):
    z = i / 1000
    dist = get_dl_and_kpc_per_asec(z=i/100)[0]
    L = (0.07 * 0.001 * 1E-26 * 4. * np.pi * dist * dist)/((1. + z) ** (1-0.7))
    lum07.append(L)


def g(x, pos):
    f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    return r'${}$'.format(f._formatSciNotation('%1.10e' % x))


def my_plot(my_directory='/home/sean/Downloads'):
    """Make some plots.

    Parameters
    ----------
    my_directory : string
        Working directory.
    """

    df = pd.read_csv(f'{my_directory}/csv.download.csv')
    # f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    # g = lambda x, pos: '${}$'.format(f._formatSciNotation('%1.10e' % x))

    plt.figure(figsize=(18, 9)).patch.set_facecolor('white')
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
    plt.errorbar(df['Extent (kpc)'], df['LDR2-to-FIRST index'],
                 xerr=df['Extent error (kpc)'],
                 yerr=df['LDR2-to-FIRST index error'],
                 marker='o', ls='none', color='k')
    plt.xlabel('Extent (kpc)', fontsize=30)
    plt.ylabel(r'$\alpha_{\mathrm{MHz-GHz}}$', fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.tight_layout()
    # plt.show()
    plt.savefig(f'{my_directory}/extent-v-index.png')
    plt.clf()

    # plot extent (kpc) against radio luminosity
    plt.errorbar(df['Extent (kpc)'], df['Luminosity with FIRST index (W/Hz)'],
                 marker='o',
                 xerr=df['Extent error (kpc)'],
                 yerr=df['Luminosity error with FIRST index (W/Hz)'],
                 ls='none', color='k')
    plt.xlabel('Extent (kpc)', fontsize=30)
    plt.ylabel(r'$L$ (W Hz$^{-1}$)', fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks([0, 1e26, 2e26, 3e26, 4e26, 5e26],
               ['0', r'$1 \times 10^{26}$', r'$2 \times 10^{26}$',
                r'$3 \times 10^{26}$', r'$4 \times 10^{26}$',
                r'$5 \times 10^{26}$'], fontsize=30)
    plt.yscale('log')
    plt.tight_layout()
    # plt.show()
    plt.savefig(f'{my_directory}/extent-v-luminosity.png')
    plt.clf()

    from matplotlib import gridspec
    # plt.margins(top=0.908,
    # bottom=0.107,
    # left=0.089,
    # right=0.947,
    # hspace=0.2,
    # wspace=0.2)
    gs = gridspec.GridSpec(2, 3, width_ratios=[4, 18, 0.4],
                           height_ratios=[18, 6])
    gs.update(wspace=0.05, hspace=0.05)
    ax = plt.subplot(gs[0, 1])
    # plot luminosity against redshift
    # s = [20*4**n for n in range(len(x))]
    s= []
    for i,j in zip(df['Compact'],df['Extent (kpc)']):
        if i:
            s.append(50)
        else:
            s.append(j**1.2)

    # s = np.array(df['Extent (kpc)'].to_list()) ** 1.2
    # print('----')
    # for a,b in zip(s, df['Extent (kpc)']):
    #     print(a,b**1.2)
    # print('----')
    c = []
    for a in df['LDR2-to-FIRST index']:
        # if a < -1:
        #     a = -1
        # print(a)
        # if a > 1:
        #     print('!!')
        #     a = 1
        c.append(a)#(a + 1) / 2)
    print(c)
    # plt.errorbar(df['redshift'], df['Luminosity with FIRST index (W/Hz)'],
    #              marker='o',
    #              #xerr=df['Redshift'] * 0.1,
    #              yerr=df['Luminosity error with FIRST index (W/Hz)'],
    #              ls='none', color='k')
    ecs=[]
    for i in df['Compact']:
        if i == False:
            ecs.append('k')
        else:
            ecs.append('gray')
    noextent = df[df['Compact']==True]
    hb=ax.scatter(df['redshift'], df['Luminosity with FIRST index (W/Hz)'],s=s,
                c=df['LDR2-to-FIRST index'], cmap='RdBu', edgecolors=ecs, linewidths=2)
    # hb=ax.scatter(df['redshift'], df['Luminosity with FIRST index (W/Hz)'],s=50,
    #             c=df['LDR2-to-FIRST index'], cmap='RdBu', edgecolors='k', linewidths=0, marker='s')
    # ax.scatter([0.4], [1e24], s=20, c=[-10], cmap='RdBu')
    ax.plot(redshift, lum, label='alpha zero', ls='dashed', color='k', lw=2)
    # plt.plot(redshift, lum07, label='alpha 7')
    ax.set_xlim(0.0, 0.7)
    ax.set_xticks([])
    # plt.xlabel(r'$z$', fontsize=30)
    # plt.ylabel(r'$L$ (W Hz$^{-1}$)', fontsize=30)
    # plt.xticks(fontsize=30)
    # plt.yticks([0, 1e26, 2e26, 3e26, 4e26, 5e26],
    #            ['0', r'$1 \times 10^{26}$', r'$2 \times 10^{26}$',
    #             r'$3 \times 10^{26}$', r'$4 \times 10^{26}$',
    #             r'$5 \times 10^{26}$'], fontsize=30)
    ax.set_yscale('log')
    ax.set_ylim(1e23, 1e27)
    ax.set_yticks([],[])
    ax.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,      # ticks along the bottom edge are off
        right=False,         # ticks along the top edge are off
        labelleft=False) # labels along the bottom edge are off

    # plt.tight_layout()

    # custom_lines = [plt.plot([0], [0], mfc='white', lw=0, marker='o', #markersize=10,
    #                        markersize=50, mec='k', mew=2, label='50 kpc'),
    #                 plt.plot([0], [0], marker='o', linestyle='None', mew=2,
    #                        markersize=100, mfc='w', mec='k', label='100 kpc')]
    # import numpy as np; np.sqrt(50**1.2)
    # https://stackoverflow.com/a/47403507/6386612
    custom_lines = [Line2D([0], [0], marker='o', linestyle='None', mew=2,
                   markersize=np.sqrt(50), mec='gray', mfc='w'),
                    Line2D([0], [0], mfc='w', mec='k',marker='o',lw=0,
                           markersize=np.sqrt(50**1.2), mew=2),
                    Line2D([0], [0], marker='o', linestyle='None', mew=2,
                           markersize=np.sqrt(500**1.2), mec='k', mfc='w')
                   ]

    # custom_labels = [f'{names[0]}s ({np.median(width[names[0]]):.0f} kpc)',
    #                  f'{names[1]}s ({np.median(width[names[1]]):.0f} kpc)']

    custom_labels = ['Unresolved','50 kpc',
                     '500 kpc']
    legend = plt.legend(custom_lines, custom_labels, ncol=3,loc='upper center',
                        bbox_to_anchor=(0.67, 1.2,-0.4,0),#loc='lower left',
                        mode='expand', numpoints=1, fontsize=30, frameon=False)
    cax = plt.subplot(gs[0, 2])
    cbar = plt.colorbar(hb, cax=cax, format='%.1f')#, vmin=-5,vmax=5)
    hb.set_clim(vmin=-1,vmax=1)
    cbar.set_label(r'$\alpha$', fontsize=30)
    # ScalarMapp
    cbar.ax.tick_params(labelsize=30)

    axx = plt.subplot(gs[1, 1])
    df_resolved = df[df['Compact'] != True]
    df_compact = df[df['Compact'] == True]
    print(df_resolved.shape)
    axx.hist(df_resolved['redshift'], histtype='step', fill=True, color=(191/256, 171/256, 37/256,0.1), edgecolor='#7a6d18', #alpha=0.1,
             lw=2, bins=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7])#, weights=weights)
    axx.hist(df_compact['redshift'], histtype='step', fill=True, color=(50/256, 118/256, 128/256,0.1), edgecolor='#054952',
             lw=2, bins=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7])#, weights=weights)
    axx.axvline(np.average(df_resolved['redshift']), color='#7a6d18', linestyle='dashed', linewidth=2)
    axx.axvline(np.average(df_compact['redshift']), color='#054952', linestyle='dashed', linewidth=2)
    axx.set_xlim(0.0, 0.7)

    axx.yaxis.tick_right()
    axx.set_xlabel(r'$z$', fontsize=30)
    axx.set_yticks([0,2,4])
    plt.setp(axx.get_xticklabels(), fontsize=30)
    plt.setp(axx.get_yticklabels(), fontsize=30)

    axy = plt.subplot(gs[0, 0])
    logbins = np.logspace(np.log10(1e23),np.log10(1e27),5)
    axy.hist(df['Luminosity with FIRST index (W/Hz)'], orientation='horizontal', histtype='step',
         fill=True, color='#faf3dd', lw=2, ec='k',bins=logbins)#, weights=weights)
    axy.axhline(np.average(df['Luminosity with FIRST index (W/Hz)']), color='k', linestyle='dashed', linewidth=2)

    axy.set_yscale('log')
    axy.set_xticks([0,3,6,9,12])
    axy.set_ylabel(r'$L$ (W Hz$^{-1}$)', fontsize=30)
    axy.set_ylim(1e23, 1e27)
    axy.set_xlim(0, 14.4)
    plt.setp(axy.get_xticklabels(), fontsize=30)
    plt.setp(axy.get_yticklabels(), fontsize=30)
    plt.show()
    # plt.savefig(f'{my_directory}/redshift-v-luminosity-colour.png')
    plt.clf()

    # plot luminosity against spectral index
    plt.errorbar(df['LDR2-to-FIRST index'],
                 df['Luminosity with FIRST index (W/Hz)'],
                 marker='o',
                 xerr=df['LDR2-to-FIRST index error'],
                 yerr=df['Luminosity error with FIRST index (W/Hz)'],
                 ls='none', color='k')
    plt.xlabel(r'$\alpha_{\mathrm{MHz-GHz}}$', fontsize=30)
    plt.ylabel(r'$L$ (W Hz$^{-1}$)', fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.yticks([0, 1e26, 2e26, 3e26, 4e26, 5e26],
               ['0', r'$1 \times 10^{26}$', r'$2 \times 10^{26}$',
                r'$3 \times 10^{26}$', r'$4 \times 10^{26}$',
                r'$5 \times 10^{26}$'], fontsize=30)
    plt.yscale('log')
    plt.tight_layout()
    # plt.show()
    plt.savefig(f'{my_directory}/index-v-luminosity.png')


def main():
    """Make some plots.
    """
    my_plot()


if __name__ == '__main__':
    main()
