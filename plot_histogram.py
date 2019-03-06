#!/usr/bin/env python2.7

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('/mnt/closet/ldr2-blazars/catalogues/ldr2-bzcat-vaneck-ned-noise.csv')

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, figsize=(12, 8))

ax1.set_xlabel('Max image flux (mJy)')
ax1.set_ylabel('Count')
ax1.hist(df['Max image flux'].dropna(), bins=np.linspace(0, 1000, 20), alpha=0.8, color='black')
print(df[df['Max image flux'] < 10000]['Max image flux'].count())

ax2.set_xlabel('Max residual (mJy)')
ax2.yaxis.set_ticks_position('none')
ax2.hist(df['Max residual'].dropna(), bins=np.linspace(0, 30, 20), alpha=0.8, color='black')
print(df[df['RMS'] < 30]['RMS'].count())

ax3.set_xlabel('RMS (mJy)')
ax3.yaxis.set_ticks_position('none')
ax3.hist(df['RMS'].dropna(), bins=np.linspace(0, 2, 20), alpha=0.8, color='black')
print(df[df['RMS'] < 0.9]['RMS'].count())
print(df['RMS'].max())

plt.tight_layout()
plt.show()

fig, axes = plt.subplots(2, 2, figsize=(10, 10))

axes[0, 0].set_xlabel('Max image flux')
axes[0, 0].set_ylabel('RMS')
axes[0, 0].loglog(df['Max image flux'].dropna(), df['RMS'].dropna(), marker='.', ls='None', alpha=0.8, color='black')

axes[0, 1].set_xlabel('Max residual')
axes[0, 1].set_ylabel('RMS')
axes[0, 1].loglog(df['Max residual'].dropna(), df['RMS'].dropna(), marker='.', ls='None', alpha=0.8, color='black')

axes[1, 0].set_xlabel('Max image flux')
axes[1, 0].set_ylabel('Max residual')
axes[1, 0].loglog(df['Max image flux'].dropna(), df['Max residual'].dropna(), marker='.', ls='None', alpha=0.8, color='black')

fig.delaxes(axes[1, 1])

plt.tight_layout()
plt.show()
