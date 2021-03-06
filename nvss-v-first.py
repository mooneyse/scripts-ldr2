#!/usr/bin/env python3

from uncertainties import ufloat
from uncertainties.umath import log, sqrt, exp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# plt.figure(figsize=(13.92, 8.60)).patch.set_facecolor('white')
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
fig, ax = plt.subplots(figsize=(13.92, 8.60))

df = pd.read_csv('/home/sean/Downloads/LDR2 and BZCAT 10_ crossmatch -'
                 ' BL Lacs(1).csv')
df = df[df['Compact']==False]  # resolved sources
print(len(df))
ax.errorbar(df['S1.4.1'], df['FINT'], label=r'Extended ($N = {}$)'.format(len(df['Compact'])),
             xerr=df['S1.4.1'] * 0.1, yerr=df['FINT'] * 0.1, color='k',
             marker='s', ls='none', mfc='#faf3dd', mew=2, mec='k', markersize=15, elinewidth=2)


df = pd.read_csv('/home/sean/Downloads/LDR2 and BZCAT 10_ crossmatch -'
                 ' BL Lacs(1).csv')

df = df[df['Compact']==True]  # compact sources

ax.errorbar(df['S1.4.1'], df['FINT'], label=r'Unresolved ($N = {}$)'.format(len(df['Compact'])),
             xerr=df['S1.4.1'] * 0.1, yerr=df['FINT'] * 0.1, color='k',
             marker='s', ls='none', mfc='#302f2c', mew=2, mec='k', markersize=15, elinewidth=2)


ax.plot([-1, 1100], [-1, 1100], ls='dashed', lw=2, color='k', marker='')

ax.set_xlim(1, 1000)
ax.set_ylim(1, 1000)
plt.setp(ax.get_xticklabels(), fontsize=30)
plt.setp(ax.get_yticklabels(), fontsize=30)
# ax.set_xticks(fontsize=30)
# ax.set_yticks(fontsize=30)

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel(r'$S_\mathrm{NVSS}$ (mJy)', fontsize=30, color='black')
ax.set_ylabel(r'$S_\mathrm{FIRST}$ (mJy)', fontsize=30, color='black')

handles, labels = ax.get_legend_handles_labels()
handles = [h[0] for h in handles]
legend = ax.legend(handles, labels, bbox_to_anchor=(0, 1.0, 1, 0), loc='lower left',  ncol=2,
                    mode='expand', numpoints=1, fontsize=30, frameon=False)
# legend = plt.legend(bbox_to_anchor=(0, 1.0, 1, 0), loc='lower left',  ncol=2,
#                     mode='expand', numpoints=1, fontsize=30, frameon=False)
plt.setp(legend.get_texts(), color='black')

#plt.show()
plt.tight_layout()
plt.savefig('/home/sean/Downloads/nvss-v-first.png')
