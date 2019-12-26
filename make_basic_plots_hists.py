#!/usr/bin/env python3

"""Make some plots.
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from matplotlib.lines import Line2D


def my_plot(my_folder='/home/sean/Downloads'):
    """Make some plots.

    Parameters
    ----------
    my_directory : string
        Working directory.
    """
    df = pd.read_csv(f'{my_folder}/csv.download.csv')
    compact = df[df['Compact'] == True]  # noqa
    extended = df[df['Compact'] != True]  # noqa
    plt.figure(figsize=(12, 12)).patch.set_facecolor('white')
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

    # luminosity against redshift
    gs = gridspec.GridSpec(2, 1, width_ratios=[1], height_ratios=[2, 1],
                           hspace=0.02)
    ax = plt.subplot(gs[0, 0])
    ax.errorbar(compact['redshift'],
                compact['Luminosity with FIRST index (W/Hz)'],
                markersize=15, mec='k', label='Unresolved', elinewidth=2,
                yerr=compact['Luminosity error with FIRST index (W/Hz)'],
                marker='s', ls='none', mfc='#302f2c', color='k', mew=2)
    ax.errorbar(extended['redshift'],
                extended['Luminosity with FIRST index (W/Hz)'],
                markersize=15, label='Extended', elinewidth=2, mec='k',
                yerr=extended['Luminosity error with FIRST index (W/Hz)'],
                marker='s', ls='none', mfc='#faf3dd', color='k', mew=2)
    ax.set_xticks([])
    ax.set_yscale('log')
    ax.set_ylabel(r'$L \,\,\,\,(\mathrm{W\,\,Hz}^{-1})$', fontsize=30)
    ax.set_xlim(0, 0.7)
    plt.setp(ax.get_yticklabels(), fontsize=30)
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] for h in handles]
    ax.legend(handles, labels, ncol=2, loc='upper center', numpoints=1,
              fontsize=30, mode='expand',
              bbox_to_anchor=(0.5, 1.2, -0.1, 0), frameon=False)

    axx = plt.subplot(gs[1, 0])
    bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    axx.hist([compact['redshift'], extended['redshift']],
             histtype='stepfilled', color=['#302f2c', '#faf3dd'],
             edgecolor='k', lw=2, bins=bins, stacked=True)
    axx.set_yticks([0, 4, 8])
    axx.set_xlabel(r'$z$', fontsize=30)
    axx.set_xlim(0, 0.7)
    plt.setp(axx.get_xticklabels(), fontsize=30)
    plt.setp(axx.get_yticklabels(), fontsize=30)
    gs.update(left=0.15)
    plt.savefig(f'{my_folder}/luminosity-against-redshift-hist.png')
    plt.clf()

    # redshift against spectral index
    gs = gridspec.GridSpec(2, 1, width_ratios=[1],
                           height_ratios=[2, 1], hspace=0.02)
    ax = plt.subplot(gs[0, 0])
    ax.errorbar(compact['LDR2-to-FIRST index'],
                compact['redshift'],
                markersize=15,
                xerr=compact['LDR2-to-FIRST index error'],
                mec='k',
                marker='s', ls='none', mfc='#302f2c', color='k', mew=2,
                label='Unresolved', elinewidth=2)
    ax.plot([0.48],
            [0.496],
            markersize=15,
            mec='k', mew=2,
            marker='>', ls='none', mfc='#302f2c', color='k')
    ax.text(0.42, 0.44, r'$\alpha = 1.95$', rotation=90, fontsize=20)
    ax.errorbar(extended['LDR2-to-FIRST index'],
                extended['redshift'],
                markersize=15,
                xerr=extended['LDR2-to-FIRST index error'],
                marker='s', ls='none', mfc='#faf3dd', color='k', mew=2,
                label='Extended', elinewidth=2, mec='k')
    ax.set_xticks([])
    ax.set_ylabel(r'$z$', fontsize=30)
    ax.set_xlim(-1, 0.5)
    plt.setp(ax.get_yticklabels(), fontsize=30)
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] for h in handles]
    ax.legend(handles, labels, ncol=2, loc='upper center', numpoints=1,
              fontsize=30, mode='expand',
              bbox_to_anchor=(0.5, 1.2, -0.1, 0), frameon=False)

    axx = plt.subplot(gs[1, 0])
    bins = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5]
    axx.hist([compact['LDR2-to-FIRST index'], extended['LDR2-to-FIRST index']],
             histtype='stepfilled', color=['#302f2c', '#faf3dd'],
             edgecolor='k', lw=2, bins=bins, stacked=True)
    axx.set_xticks([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5])
    axx.set_yticks([0, 4, 8])
    axx.set_xlabel(r'$\alpha$', fontsize=30)
    axx.set_xlim(-1, 0.5)
    plt.setp(axx.get_xticklabels(), fontsize=30)
    plt.setp(axx.get_yticklabels(), fontsize=30)
    gs.update(left=0.15)
    plt.savefig(f'{my_folder}/redshift-against-spectral-index-hist.png')
    plt.clf()

    # spectral index against luminosity
    gs = gridspec.GridSpec(2, 1, width_ratios=[1],
                           height_ratios=[2, 1], hspace=0.02)
    ax = plt.subplot(gs[0, 0])
    ax.errorbar(compact['Luminosity with FIRST index (W/Hz)'],
                compact['LDR2-to-FIRST index'],
                markersize=15,
                xerr=compact['Luminosity error with FIRST index (W/Hz)'],
                yerr=compact['LDR2-to-FIRST index error'],
                mec='k',
                marker='s', ls='none', mfc='#302f2c', color='k', mew=2,
                label='Unresolved', elinewidth=2)
    ax.plot([1.31e24],
            [0.46],
            markersize=15,
            mec='k', mew=2,
            marker='^', ls='none', mfc='#302f2c', color='k')
    ax.text(0.9e24, 0.37, r'$\alpha = 1.95$', rotation=00, fontsize=20)
    ax.errorbar(extended['Luminosity with FIRST index (W/Hz)'],
                extended['LDR2-to-FIRST index'],
                markersize=15,
                xerr=extended['Luminosity error with FIRST index (W/Hz)'],
                yerr=extended['LDR2-to-FIRST index error'],
                marker='s', ls='none', mfc='#faf3dd', color='k', mew=2,
                label='Extended', elinewidth=2, mec='k')
    bins = np.logspace(np.log10(1e23), np.log10(1e27), 9)
    ax.set_ylabel(r'$\alpha$', fontsize=30)
    ax.set_xlim(bins[1], bins[-1])
    ax.set_xscale('log')
    ax.set_xticks([])
    ax.set_ylim(-1, 0.5)
    ax.tick_params(axis='x',
                   which='both',
                   bottom=False,
                   top=False,
                   labelbottom=False)
    plt.setp(ax.get_yticklabels(), fontsize=30)
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] for h in handles]
    ax.legend(handles, labels, ncol=2, loc='upper center', numpoints=1,
              fontsize=30, mode='expand',
              bbox_to_anchor=(0.5, 1.2, -0.1, 0), frameon=False)

    axx = plt.subplot(gs[1, 0])
    axx.hist([compact['Luminosity with FIRST index (W/Hz)'],
             extended['Luminosity with FIRST index (W/Hz)']],
             histtype='stepfilled', color=['#302f2c', '#faf3dd'],
             edgecolor='k', lw=2, stacked=True, bins=bins)
    axx.set_yticks([0, 4, 8])
    axx.set_xlabel(r'$L \,\,\,\,(\mathrm{W\,\,Hz}^{-1})$', fontsize=30)
    axx.set_xlim(bins[1], bins[-1])
    axx.set_xscale('log')
    plt.setp(axx.get_xticklabels(), fontsize=30)
    plt.setp(axx.get_yticklabels(), fontsize=30)
    gs.update(left=0.15)
    plt.savefig(f'{my_folder}/spectral-index-against-luminosity-hist.png')
    plt.close()

    '''---------------------------------------------------------------------'''
    # # spectral index against luminosity
    # gs = gridspec.GridSpec(2, 1, width_ratios=[1],
    #                        height_ratios=[2, 1], hspace=0.02)
    # ax = plt.subplot(gs[0, 0])
    # # ax.errorbar(compact['Core dominance'],
    # #             compact['Extent (kpc)'],
    # #             markersize=15,
    # #             xerr=compact['Core dominance error'],
    # #             yerr=compact['Extent error (kpc)'],
    # #             mec='k',
    # #             marker='s', ls='none', mfc='#302f2c', color='k', mew=2,
    # #             label='Unresolved', elinewidth=2)
    # # ax.plot([1.31e24],
    # #         [0.46],
    # #         markersize=15,
    # #         mec='k', mew=2,
    # #         marker='^', ls='none', mfc='#302f2c', color='k')
    # # ax.text(0.9e24, 0.37, r'$\alpha = 1.95$', rotation=00, fontsize=20)
    # ax.errorbar(extended['Core dominance'],
    #             extended['Extent (kpc)'],
    #             markersize=15,
    #             xerr=extended['Core dominance error'],
    #             yerr=extended['Extent error (kpc)'],
    #             marker='s', ls='none', mfc='#faf3dd', color='k', mew=2,
    #             label='Extended', elinewidth=2, mec='k')
    # bins = np.logspace(np.log10(1e23), np.log10(1e27), 9)
    # ax.set_ylabel(r'$R$', fontsize=30)
    # # ax.set_xlim(bins[1], bins[-1])
    # # ax.set_xscale('log')
    # ax.set_xticks([])
    # # ax.set_ylim(-1, 0.5)
    # ax.tick_params(axis='x',
    #                which='both',
    #                bottom=False,
    #                top=False,
    #                labelbottom=False)
    # plt.setp(ax.get_yticklabels(), fontsize=30)
    # handles, labels = ax.get_legend_handles_labels()
    # handles = [h[0] for h in handles]
    # ax.legend(handles, labels, ncol=2, loc='upper center', numpoints=1,
    #           fontsize=30, mode='expand',
    #           bbox_to_anchor=(0.5, 1.2, -0.1, 0), frameon=False)
    #
    # axx = plt.subplot(gs[1, 0])
    # axx.hist(#[compact['Core dominance'],
    #          extended['Core dominance'],#],
    #          histtype='stepfilled', color=['#302f2c', '#faf3dd'],
    #          edgecolor='k', lw=2, stacked=True, bins=bins)
    # axx.set_yticks([0, 4, 8])
    # axx.set_xlabel(r'$L \,\,\,\,(\mathrm{W\,\,Hz}^{-1})$', fontsize=30)
    # axx.set_xlim(bins[1], bins[-1])
    # axx.set_xscale('log')
    # plt.setp(axx.get_xticklabels(), fontsize=30)
    # plt.setp(axx.get_yticklabels(), fontsize=30)
    # gs.update(left=0.15)
    # plt.savefig(f'{my_folder}/core-dominance-against-extent-hist.png')
    # plt.close()





    '''---------------------------------------------------------------------'''

    plt.figure(figsize=(20, 12)).patch.set_facecolor('white')
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
    gs = gridspec.GridSpec(2, 3, width_ratios=[4, 18, 0.4],
                           height_ratios=[18, 6])
    gs.update(wspace=0.05, hspace=0.05)
    ax = plt.subplot(gs[0, 1])
    extended.sort_values('Luminosity with FIRST index (W/Hz)', inplace=True,
                         ascending=False)
    x = np.log10(extended['Luminosity with FIRST index (W/Hz)'])
    s = (((60 - 10) / (np.max(x) - np.min(x))) * (x - np.max(x)) + 60) ** 2
    hb = ax.scatter(extended['Extent (kpc)'].astype('float'),
                    extended['Core dominance'].astype('float'),
                    s=s, cmap='plasma', linewidths=2, edgecolors='k',
                    c=extended['LDR2-to-FIRST index'])
    # ax.set_yscale('log')
    ax.set_xlim(0, 700)
    ax.set_ylim(0, 1)#1e-2, 1e1)
    ax.set_xticks([])
    ax.tick_params(axis='y',
                   which='both',
                   left=False,
                   right=False,
                   labelleft=False)
    # ax.plot()
    # ax.text(0.42, 0.44, r'$R = 1726$', rotation=0, fontsize=20)

    custom_lines = [Line2D([0], [0], mfc='w', mec='k', marker='o', lw=0,
                           markersize=np.sqrt(s.iloc[0]), mew=2),
                    Line2D([0], [0], marker='o', linestyle='None', mew=2,
                           markersize=np.sqrt(s.iloc[-1]), mec='k', mfc='w')]
    labels = [f'{extended["Luminosity with FIRST index (W/Hz)"].iloc[0]}',
              f'{extended["Luminosity with FIRST index (W/Hz)"].iloc[-1]}']
    labels = [r'$1.9 \times 10^{24}$ W Hz', r'$4.5 \times 10^{26}$ W Hz']
    plt.legend(custom_lines, labels, ncol=2, loc='upper center',
               bbox_to_anchor=(0.55, 1.2, -0.15, 0),
               mode='expand', numpoints=1, fontsize=30, frameon=False)

    cax = plt.subplot(gs[0, 2])
    cbar = plt.colorbar(hb, cax=cax, format='%.1f')
    hb.set_clim(vmin=-1, vmax=0)
    cbar.set_label(r'$\alpha$', fontsize=30)
    cbar.ax.tick_params(labelsize=30)

    axx = plt.subplot(gs[1, 1])
    axx.hist(extended['Extent (kpc)'].astype('float'), histtype='step',
             fill=True, color='#faf3dd',
             edgecolor='k',
             lw=2, bins=[0, 100, 200, 300, 400, 500, 600, 700])
    axx.axvline(np.average(extended['Extent (kpc)']), color='k',
                linewidth=2, linestyle='dashed')
    axx.set_xlim(0, 700)
    axx.yaxis.tick_right()
    axx.set_xlabel(r'$D$ (kpc)', fontsize=30)
    axx.set_yticks([0, 3, 6])
    plt.setp(axx.get_xticklabels(), fontsize=30)
    plt.setp(axx.get_yticklabels(), fontsize=30)

    axy = plt.subplot(gs[0, 0])
    logbins = np.logspace(np.log10(1e-2), np.log10(1e4), 13)
    axy.hist(extended['Core dominance'].astype('float'),
             orientation='horizontal',
             histtype='step', fill=True, color='#faf3dd', lw=2, ec='k',
             bins=logbins)
    axy.axhline(np.median(extended['Core dominance'].astype('float')),
                color='k',
                linestyle='dashed', linewidth=2)
    axy.set_yscale('log')
    axy.set_xticks([0, 3, 6, 9])
    axy.set_ylabel(r'$R$', fontsize=30)
    axy.set_ylim(1e-2, 1e1)
    plt.setp(axy.get_xticklabels(), fontsize=30)
    plt.setp(axy.get_yticklabels(), fontsize=30)
    # gs.update(bottom=0)
    # plt.show()
    print(len(extended['Core dominance']))
    plt.savefig(f'{my_folder}/core-dominance-against-extent-hist.png')


def main():
    """Make some plots.
    """
    my_plot()


if __name__ == '__main__':
    main()
