#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import os


def luminosity(s, d, z, alpha):
    return (s * 4 * np.pi * (d ** 2)) / ((1 + z) ** (1 + alpha))


df = pd.read_csv('/home/sean/Downloads/'
                 'LDR2 and BZCAT 10_ crossmatch - BL Lacs.csv')
L_total_cpt, L_total_ext, L_core, L_ext = [], [], [], []
for s_total, s_core, s_ext, d, z, alpha, cpt in zip(df['Total_flux'],
                                                    df['Core flux (mJy)'],
                                                    df['Diffuse MHz'],
                                                    df['Luminosity distance'],
                                                    df['redshift'],
                                                    df['LDR2-to-NVSS index'],
                                                    df['Compact']):
    if cpt:
        L_total_cpt.append(luminosity(s=s_total * 1e-3 * 1e-26, d=d, z=z,
                                      alpha=alpha))
    else:
        L_total_ext.append(luminosity(s=s_total * 1e-3 * 1e-26, d=d, z=z,
                                      alpha=alpha))
        L_core.append(luminosity(s=s_core * 1e-3 * 1e-26, d=d, z=z, alpha=0))
        L_ext.append(luminosity(s=(s_total - s_core) * 1e-3 * 1e-26, d=d, z=z,
                                alpha=-0.8))

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
fig = plt.figure(figsize=(14, 20))
bins = np.logspace(np.log10(1e23), np.log10(1e27), 9)
color = '#FAF3DD'
lw = 3
ec = 'k'
histtype = 'step'
fontsize = 30

ax1 = plt.subplot(311)
plt.hist([L_total_ext, L_total_cpt],
         bins=bins,
         density=False,
         cumulative=False,
         color=[color, '#302F2C'],
         lw=lw,
         # fill=True,
         histtype='stepfilled',
         stacked=True,
         label=['Extended', 'Unresolved'],
         ec=ec)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.tick_params(axis='x', direction='in', top=True, bottom=True, which='both')
plt.yticks(np.linspace(2, 12, 6), fontsize=fontsize)
plt.xscale('log')
plt.xlim(bins[0], bins[-1])
plt.ylim(0, 13)
plt.text(2.5e26, 8, r'$L_\mathrm{total}$', fontsize=fontsize)
plt.ylabel(r'$N$', fontsize=30)
legend = plt.legend(ncol=2, fontsize=fontsize, frameon=False, mode='expand',
                    bbox_to_anchor=(0,0.2,1,1))

# )#bbox_to_anchor=(0, 0, 0, 0), loc='upper left',
                    #mode='expand', numpoints=1, fontsize=30, frameon=False,
                    #ncol=2)
plt.setp(legend.get_texts(), color='black')

ax2 = plt.subplot(312, sharex=ax1, sharey=ax1)
plt.hist(L_core,
         bins=bins,
         density=False,
         cumulative=False,
         color=color,
         lw=lw,
         fill=True,
         histtype=histtype,
         ec=ec)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.tick_params(axis='x', direction='in', top=True, bottom=True, which='both')
plt.yticks(fontsize=fontsize)
plt.text(2.5e26, 8, r'$L_\mathrm{core}$', fontsize=fontsize)
plt.ylabel(r'$N$', fontsize=30)

ax3 = plt.subplot(313, sharex=ax1, sharey=ax1)
plt.hist(L_ext,
         bins=bins,
         density=False,
         cumulative=False,
         color=color,
         lw=lw,
         fill=True,
         histtype=histtype,
         ec=ec)
plt.tick_params(axis='x', direction='in', top=True, bottom=True, which='both',
                pad=15)
plt.yticks(fontsize=fontsize)
plt.text(2.5e26, 8, r'$L_\mathrm{ext}$', fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.ylabel(r'$N$', fontsize=30)

plt.xlabel(r'$L_{144}$ (W Hz$^{-1}$)', fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.subplots_adjust(wspace=0, hspace=0)
save_name = '/home/sean/Downloads/luminosity-hist.png'
plt.savefig(save_name)
os.system(f'convert {save_name} -trim {save_name}')  # removes whitespace
print(f'gpicview {save_name}')
