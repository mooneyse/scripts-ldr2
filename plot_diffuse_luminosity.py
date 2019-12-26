#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from uncertainties import ufloat


def values_with_errors(value, error):
    '''Return a list with uncertainties on the floats.
    '''
    combined = []
    for v, e in zip(value, error):
        combined.append(ufloat(v, e))
    return np.array(combined)


def l_nu(s, alpha, z, d):
    '''Calculate the radio luminoisty in watts per hertz.
    '''
    s = 1e-26 * 1e-3 * s  # millijansky to watts
    return (s * 4 * np.pi * ( d ** 2)) / ((1 + z) ** ( 1 + alpha))


# get the data
df = pd.read_csv('/home/sean/Downloads/csv.download.csv')  # read data
df = df[df['Compact'] == False]  # filter out unresolved sources

l = df['Luminosity with NVSS index (W/Hz)']  # luminosity
l_error = df['Luminosity error with NVSS index (W/Hz)']  # luminosity
l_total = values_with_errors(value=l, error=l_error)

s = df['Total_flux']  # total 144 mhz flux density
s_error = df['E_Total_flux']  # total 144 mhz flux density uncertainty
s_total = values_with_errors(value=s, error=s_error)

s = df['Core flux (mJy)']
s_error =df['Core flux error (mJy)']
s_core = values_with_errors(value=s, error=s_error)

s_diffuse = s_total - s_core

alpha_core = ufloat(0, 0.2)
alpha_diffuse = ufloat(-0.8, 0.2)  # spectral index of the diffuse component

d = df['Luminosity distance']  # luminosity distance
z = df['redshift']  # redshift with no uncertainty

l_core, l_diffuse = [], []
for core, diffuse, redshift, l_distance in zip(s_core, s_diffuse, z, d):
    l_core.append(l_nu(s=core, alpha=alpha_core, z=redshift, d=l_distance))
    l_diffuse.append(l_nu(s=diffuse, alpha=alpha_diffuse, z=redshift,
                          d=l_distance))

l_core = np.array(l_core)
l_diffuse = np.array(l_diffuse)

# print the results
print('Source, Core flux, Core luminosity, Diffuse flux, Diffuse luminosity')
for s, fc, lc, fd, ld in zip (df['Name'], s_core, l_core, s_diffuse,
                              l_diffuse):
    print(s, fc, lc, fd, ld)

print(f'Average L_tot:  {np.mean(l_total):.2E}')
print(f'Average L_core: {np.mean(l_diffuse):.2E}')
print(f'Average L_ext: {np.mean(l_core):.2E}')

# make the plot
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

low = 1e23
high = 1e27
bins = np.logspace(np.log10(low), np.log10(high), (np.log10(high) -
                                                   np.log10(low)) * 2 + 1)
l_total = [l_total[i].n for i, _ in enumerate(l_total)]
l_core = [l_core[i].n for i, _ in enumerate(l_core)]
l_diffuse = [l_diffuse[i].n for i, _ in enumerate(l_diffuse)]

plt.hist(l_total, bins=bins, histtype='step', lw=5, fill=True, label='Total',
         ec=(52/255, 84/255, 209/255), fc=(52/255, 84/255, 209/255, 0.1))

plt.hist(l_core, bins=bins, histtype='step', label='Core', lw=5, fill=True,
         ec=(209/255, 52/255, 91/255), fc=(209/255, 52/255, 91/255, 0.1),
         ls='dashed')

plt.hist(l_diffuse, bins=bins, histtype='step', label='Extended', lw=5,
         ec=(247/255, 203/255, 21/255), fc=(247/255, 203/255, 21/255, 0.1),
         ls='dotted', fill=True)

plt.xscale('log')
plt.xlim(low, high)
plt.xlabel(r'$L_{144}$', fontsize=30)
plt.ylabel(r'$N$', fontsize=30)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

legend = plt.legend(bbox_to_anchor=(0, 1, 1, 0), loc='lower left',  ncol=3,
                    mode='expand', numpoints=1, fontsize=30, frameon=False)
plt.setp(legend.get_texts(), color='black')
plt.tight_layout()
plt.savefig('/home/sean/Downloads/extended-luminosity.png')

# bring in uncertainties