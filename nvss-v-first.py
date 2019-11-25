#!/usr/bin/env python3

from uncertainties import ufloat
from uncertainties.umath import log, sqrt, exp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

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

df = pd.read_csv('/home/sean/Downloads/LDR2 and BZCAT 10_ crossmatch -'
                 ' Sheet7.csv')
df = df[df['Resolved?']==True]  # compact sources

plt.errorbar(df['S1.4.1'], df['FINT'], label=r'Unresolved ($N = {}$)'.format(len(df['Resolved?'])),
             xerr=df['S1.4.1'] * 0.1, yerr=df['FINT'] * 0.1, color='k',
             marker='s', ls='none', mfc='k', mew=2, mec='k', markersize=10)

df = pd.read_csv('/home/sean/Downloads/LDR2 and BZCAT 10_ crossmatch -'         
                 ' Sheet7.csv')                                                 
df = df[df['Resolved?']==False]  # resolved sources                                               

plt.errorbar(df['S1.4.1'], df['FINT'], label=r'Extended ($N = {}$)'.format(len(df['Resolved?'])),
             xerr=df['S1.4.1'] * 0.1, yerr=df['FINT'] * 0.1, color='k',
             marker='s', ls='none', mfc='w', mew=2, mec='k', markersize=10)  

plt.plot([-1, 1100], [-1, 1100], ls='dashed', lw=2, color='k', marker='')

plt.xlim(1, 1000)
plt.ylim(1, 1000)

plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

plt.xscale('log')
plt.yscale('log')

plt.xlabel('NVSS flux density (mJy)', fontsize=30, color='black')
plt.ylabel('FIRST flux density (mJy)', fontsize=30, color='black')

legend = plt.legend(bbox_to_anchor=(0, 1.0, 1, 0), loc='lower left',  ncol=2,
                    mode='expand', numpoints=1, fontsize=30, frameon=False)
plt.setp(legend.get_texts(), color='black')

#plt.show()
plt.tight_layout()
plt.savefig('/home/sean/Downloads/nvss-v-first.png')