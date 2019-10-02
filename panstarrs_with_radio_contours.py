#!/usr/bin/env python3

"""Plot Pan-STARRS data.
"""

import matplotlib as mpl
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from ds9norm import DS9Normalize
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from bl_lac_postage_stamps import get_kpc_per_asec

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '06 September 2019'


def optical(sigma=4, my_directory='/data5/sean/ldr2'):
    """Plot Pan-STARRS data.

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
    pix = 0.25  # arcseconds per pixel

    for source_name, ra, dec, mosaic, rms, z in zip(df['Source name'],
                                                    df['RA (J2000.0)'],
                                                    df['Dec (J2000.0)'],
                                                    df['Mosaic_ID'],
                                                    df['Isl_rms'],
                                                    df['Redshift']):

        sky_position = SkyCoord(ra, dec, unit='deg')
        if source_name == '5BZBJ1202+4444':
            size = [3, 3] * u.arcmin
            p = 54
        elif source_name == '5BZBJ1419+5423':
            size = [4, 4] * u.arcmin
            p = 72
        else:
            size = [2, 2] * u.arcmin
            p = 36

        p_hdu = fits.open(f'{my_directory}/panstarrs/{source_name}.i.fits')[0]
        p_wcs = WCS(p_hdu.header)
        p_cutout = Cutout2D(p_hdu.data, sky_position, size=size, wcs=p_wcs)

        l_hdu = fits.open(f'{my_directory}/mosaics/{mosaic}-mosaic.fits')[0]
        l_wcs = WCS(l_hdu.header, naxis=2)
        l_cutout = Cutout2D(np.squeeze(l_hdu.data), sky_position, size=size,
                            wcs=l_wcs)

        levels = [level * rms / 1000 for level in [4, 8, 16, 32]]

        ax = plt.subplot(projection=p_cutout.wcs)
        im = ax.imshow(p_cutout.data, vmin=0, vmax=8000, cmap='Greys',
                       origin='lower', norm=DS9Normalize(stretch='arcsinh'))

        cbar = plt.colorbar(im)
        cbar.set_label('Excess counts', size=20)
        cbar.ax.tick_params(labelsize=20)

        ax.contour(l_cutout.data, transform=ax.get_transform(l_cutout.wcs),
                   levels=levels, origin='lower', colors=colors)

        plt.xlabel('Right ascension', fontsize=20, color='black')
        plt.ylabel('Declination', fontsize=20, color='black')
        ax.tick_params(axis='both', which='major', labelsize=20)
        plt.minorticks_on()
        ax.tick_params(which='minor', length=0)

        kpc_per_asec = get_kpc_per_asec(z=z)
        sbar = sbar_asec / pix  # length of scalebar in pixels
        kpc_per_pixel = kpc_per_asec * pix
        s = p_cutout.data.shape[1]  # plot scalebar
        plt.plot([p, p + sbar], [s - p, s - p], marker='None', lw=2, color='b')
        plt.text(p,  s - p + 10 * (p / 36), f'{sbar_asec:.0f}" = '
                 f'{kpc_per_pixel * sbar:.0f} kpc', fontsize=20, color='b')

        plt.savefig(f'{my_directory}/images/panstarrs-{source_name}.png')
        plt.clf()


def main():
    """Plot Pan-STARRS data.
    """
    optical()


if __name__ == '__main__':
    main()
