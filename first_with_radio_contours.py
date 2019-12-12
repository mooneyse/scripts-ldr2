#!/usr/bin/env python3

"""Plot FIRST data.
"""

import matplotlib as mpl
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from ds9norm import DS9Normalize
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from bl_lac_postage_stamps import get_kpc_per_asec

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '06 September 2019'


def ghz(sigma=4, my_directory='/data5/sean/ldr2'):
    """Plot FIRST data.

    Parameters
    ----------
    sigma : float or int
        The threshold of the significance to set the mask, as a factor of the
        local RMS. The default is 4.
    my_directory : string
        Working directory.
    """
    df = pd.read_csv(f'{my_directory}/catalogues/final.csv')

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
    colors = ['#118ab2', '#06d6a0', '#ffd166', '#ef476f']
    sbar_asec = 30  # desired length of scalebar in arcseconds
    pix = 1.8  # arcseconds per pixel

    for source_name, ra, dec, mosaic, rms, z, f_rms, pf in zip(
        df['Name'], df['BZCAT RA'], df['BZCAT Dec'],
        df['Mosaic_ID'], df['Isl_rms'], df['redshift'],
            df['RMS'], df['Peak_flux']):
        source_name = '5BZB' + source_name
        sky_position = SkyCoord(ra, dec, unit='deg')
        if source_name == '5BZBJ1202+4444' or source_name == '5BZBJ1325+4115':
            size = [3, 3] * u.arcmin
            x = 9 * (1.5 / pix)  # to get the scalebar in the same place
        elif (source_name == '5BZBJ1419+5423' or
              source_name == '5BZBJ0945+5757'):
            size = [4, 4] * u.arcmin
            x = 12 * (1.5 / pix)
        else:
            size = [2, 2] * u.arcmin
            x = 6 * (1.5 / pix)

        f_hdu = fits.open(f'{my_directory}/first/{source_name}.fits')[0]
        f_wcs = WCS(f_hdu.header, naxis=2)
        print(f'{source_name}')
        f_cutout = Cutout2D(np.squeeze(f_hdu.data), sky_position, size=size,
                            wcs=f_wcs)

        l_hdu = fits.open(f'{my_directory}/mosaics/{mosaic}-mosaic.fits')[0]
        l_wcs = WCS(l_hdu.header, naxis=2)
        l_cutout = Cutout2D(np.squeeze(l_hdu.data), sky_position, size=size,
                            wcs=l_wcs)

        if 4 * rms < (pf / 50):
            # see 2.2 of https://arxiv.org/pdf/1907.03726.pdf
            levels = [level * (pf / 50) / 1000 for level in [1, 2, 4, 8]]
        else:
            levels = [level * rms / 1000 for level in [4, 8, 16, 32]]
        f_level = f_rms * 4 / 1000

        ax = plt.subplot(projection=f_cutout.wcs)
        im = ax.imshow(f_cutout.data, vmin=0, vmax=np.max(f_cutout.data),
                       cmap='cubehelix_r', origin='lower',
                       interpolation='gaussian',
                       norm=DS9Normalize(stretch='arcsinh'))
        ax.contour(f_cutout.data, levels=[f_level], origin='lower',
                   colors=['k'])

        cbar = plt.colorbar(im)
        cbar.set_label(r'Jy beam$^{-1}$', size=20)
        cbar.ax.tick_params(labelsize=20)

        ax.contour(l_cutout.data, transform=ax.get_transform(l_cutout.wcs),
                   levels=levels, origin='lower', colors=colors)

        beam = Circle((6, 6), radius=1.5, linestyle='dashed', lw=2, fc='none',
                      edgecolor='blue')  # beam = 5.4" diameter, 1 pixel = 1.8"
        ax.add_patch(beam)

        plt.xlabel('Right ascension', fontsize=20, color='black')
        plt.ylabel('Declination', fontsize=20, color='black')
        ax.tick_params(axis='both', which='major', labelsize=20)
        plt.minorticks_on()
        ax.tick_params(which='minor', length=0)
        plt.xlim(0, f_cutout.data.shape[0])
        plt.ylim(0, f_cutout.data.shape[1])

        kpc_per_asec = get_kpc_per_asec(z=z)
        sbar = sbar_asec / pix  # length of scalebar in pixels
        kpc_per_pixel = kpc_per_asec * pix
        s = f_cutout.data.shape[1]
        plt.plot([x, x + sbar], [s - x, s - x], marker='None', lw=2, color='b')
        plt.text(x, s - (4.167 * x / 5), f'{sbar_asec:.0f}" = '
                 f'{kpc_per_pixel * sbar:.0f} kpc', fontsize=20, color='b')

        plt.savefig(f'{my_directory}/images/first-{source_name}.png')
        plt.clf()


def main():
    """Plot FIRST data.
    """
    ghz()


if __name__ == '__main__':
    main()
