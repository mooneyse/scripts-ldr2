#!/usr/bin/env python3

"""Plot postage stamp images of LDR2 BL Lacs."""

import os
import operator
import numpy as np
import pandas as pd
from math import sqrt, exp, sin
from skimage.measure import label
import matplotlib as mpl
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from ds9norm import DS9Normalize
import smallestenclosingcircle

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '06 September 2019'


def nearest_to_centre(my_arr, percent):
    """Given a two dimensional array, return the value of the pixel nearest to
    the centre that is non-zero.

    Paramters
    ---------
    my_arr : NumPy array
        Array to use.
    percent : float
        Fraction of the array to consider the centre.

    Returns
    ------
    Indices of the islands.
    """
    if np.max(my_arr) == 1:
        return [1]
    i1 = int(round(my_arr.shape[0] * (0.5 - percent / 2), 0))
    i2 = int(round(my_arr.shape[0] * (0.5 + percent / 2), 0))
    islands_in_the_sun = list(np.unique(my_arr[i1:i2, i1:i2]))
    if 0 in islands_in_the_sun:
        islands_in_the_sun.remove(0)

    if len(islands_in_the_sun) >= 1:
        return islands_in_the_sun
    else:  # no islands so find the nearest to the centre
        R = int(round(my_arr.shape[0] / 2, 0))
        C = int(round(my_arr.shape[0] / 2, 0))

        dist = {}
        for r in range(my_arr.shape[0]):
            for c in range(my_arr.shape[1]):
                if my_arr[r, c] != 0:
                    dist[my_arr[r, c]] = np.sqrt((R - r) ** 2 + (C - c) ** 2)
        return [min(dist.items(), key=operator.itemgetter(1))[0]]


def get_kpc_per_asec(z, H0=70, WM=0.26, WV=0.74):
    """Ned Wright's cosmology calculator. See
    http://www.astro.ucla.edu/~wright/CosmoCalc.html for the online version.

    Parameters
    ----------
    z : float
        Redshift.
    H0 : float, optional
        The Hubble constant. The default is 70. What are the units?
    WM : float, optional
        The omega matter value. The default is 0.26.
    WV : float, optional
        The omega vacuum parameter. The default is 0.74.
    n : float, optional
        The number of points to use in the integration.

    Returns
    -------
    float
        The kpc per arcsecond at the given redshift.
    """
    WR = 0  # omega(radiation)
    WK = 0  # omega curvaturve = 1 - omega(total)
    c = 299792.458  # velocity of light in km / s
    DCMR = 0  # comoving radial distance in units of c / H0
    h = H0 / 100
    WR = 4.165E-5 / (h * h)  # with 3 massless neutrino species, T0 = 2.72528
    WK = 1 - WM - WR - WV
    az = 1 / (1 + z)
    n = 1000  # number of points in integrals
    for i in range(n):
        a = az + (1 - az) * (i + 0.5) / n
        adot = sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
        DCMR = DCMR + 1 / (a * adot)
    DCMR = (1 - az) * DCMR / n
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
        ratio = 1 + y / 6 + y * y / 120
    return ((c / H0) * (az * (ratio * DCMR))) / 206.264806


def loop_through_sources(sigma=4, my_directory='/data5/sean/ldr2'):
    """Plot postage stamp images of LDR2 BL Lacs.

    Parameters
    ----------
    sigma : float or integer
        The threshold of the significance to set the mask, as a factor of the
        local RMS. The default is 4.
    my_directory : string
        Working directory.

    Returns
    -------
    string
        The name of the CSV containing the results.
    """
    results_csv = f'{my_directory}/images/ldr2.csv'
    try:
        os.remove(results_csv)
    except OSError:
        pass
    df = pd.read_csv(f'{my_directory}/images/LDR2 and BZCAT 10_ crossmatch -'
                     ' Sheet6.csv')
    result_header = ('Name,RA,Dec,RMS (uJy),Redshift,Width ("),Width (kpc),'
                     'Peak flux (mJy)\n')

    with open(results_csv, 'a') as f:
        f.write(result_header)

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
    dummy = 123456
    sbar_asec = 30  # desired length of scalebar in arcseconds
    pix = 1.5  # arcseconds per pixel

    for source_name, ra, dec, mosaic, rms, z, pf in zip(df['BZCAT name'],
                                                        df['BZCAT RA'],
                                                        df['BZCAT Dec'],
                                                        df['Mosaic_ID'],
                                                        df['Isl_rms'],
                                                        df['redshift'],
                                                        df['Peak_flux']):
        source_name = source_name.replace(' ', '')
        threshold = sigma * rms / 1000   # jansky
        thresh_ans = f'{sigma}sigma'
        if threshold < (pf / 50) / 1000:
            # see 2.2 of https://arxiv.org/pdf/1907.03726.pdf
            threshold = (pf / 50) / 1000
            thresh_ans = '1/50 S_peak'
        hdu = fits.open(f'{my_directory}/mosaics/{mosaic}-mosaic.fits')[0]
        wcs = WCS(hdu.header, naxis=2)
        sky_position = SkyCoord(ra, dec, unit='deg')
        if source_name == '5BZBJ1202+4444' or source_name == '5BZB J1202+4444':
            size = [3, 3] * u.arcmin
            p = 9
        elif source_name == '5BZBJ1419+5423':
            size = [4, 4] * u.arcmin
            p = 12
        else:
            size = [2, 2] * u.arcmin
            p = 6

        cutout = Cutout2D(np.squeeze(hdu.data), sky_position, size=size,
                          wcs=wcs)
        d = cutout.data
        copy_d = np.copy(d)
        another_copy_d = np.copy(d)
        d[d < threshold] = 0
        d[d >= threshold] = 1
        rows, cols = d.shape

        d = label(d)  # label islands of emission
        source_islands = nearest_to_centre(d, percent=0.1)
        for source_island in source_islands:
            d[d == source_island] = dummy

        d[d != dummy] = 0
        copy_d[d != dummy] = 0
        set_to_nil = []  # identify values we can set to zero for being inside
        for r in range(rows):  # set to 0 if surrounded by non nans
            for c in range(cols):
                try:
                    if (d[r - 1, c - 1] != 0 and d[r - 1, c] != 0 and
                        d[r - 1, c + 1] != 0 and d[r, c - 1] != 0 and
                        d[r, c + 1] != 0 and d[r + 1, c - 1] != 0 and
                            d[r + 1, c] != 0 and d[r + 1, c + 1] != 0):
                        set_to_nil.append((r, c))
                except IndexError:
                    print(f'Index error for {source_name}.')
                    continue

        for r, c in set_to_nil:
            d[r, c] = 0  # needs separate loop to avoid checkered pattern

        # d is an outline of the source (one) and everything else is zero
        # copy_d is the source with flux values and everything else is zero
        # another_copy_d has flux values throughout
        good_cells = []
        for r in range(rows):
            for c in range(cols):
                if d[r, c] != 0:
                    good_cells.append([r, c])

        x, y, r = smallestenclosingcircle.make_circle(good_cells)

        ax = plt.subplot(projection=cutout.wcs)
        plt.xlabel('Right ascension', fontsize=20, color='black')
        plt.ylabel('Declination', fontsize=20, color='black')
        ax.tick_params(axis='both', which='major', labelsize=20)
        plt.imshow(another_copy_d, vmin=0, vmax=np.nanmax(another_copy_d),
                   origin='lower', norm=DS9Normalize(stretch='arcsinh'),
                   cmap='plasma_r')  # interpolation='gaussian'
        beam = Circle((6, 6), radius=2, linestyle='dashed', lw=2, fc='none',
                      edgecolor='blue')  # radius=2 pixels -> 3" -> diameter=6"
        diffuse = Circle((y + 0.5, x + 0.5), radius=r, fc='none',
                         edgecolor='k', lw=2)
        ax.add_patch(beam)
        ax.add_patch(diffuse)

        kpc_per_asec = get_kpc_per_asec(z=z)
        sbar = sbar_asec / pix  # length of scalebar in pixels
        kpc_per_pixel = kpc_per_asec * pix
        s = cutout.data.shape[1]  # plot scalebar
        plt.plot([p, p + sbar], [s - p, s - p], marker='None', lw=2, color='b')
        plt.text(p, s - (5 * p / 6), f'{sbar_asec:.0f}" = '
                 f'{kpc_per_pixel * sbar:.0f} kpc', fontsize=20, color='b')

        cbar = plt.colorbar()
        cbar.set_label(r'Jy beam$^{-1}$', size=20)
        cbar.ax.tick_params(labelsize=20)
        plt.minorticks_on()
        plt.tick_params(which='minor', length=0)
        plt.contour(another_copy_d, levels=[threshold], origin='lower',
                    colors='w')
        plt.contour(another_copy_d - copy_d, levels=[threshold], colors='grey',
                    origin='lower')
        plt.savefig(f'{my_directory}/images/{source_name}.png')
        plt.clf()

        width = r * kpc_per_pixel * 2  # radius to diameter

        result = (f'{source_name},{ra},{dec},{rms * 1e3},{z},'
                  f'{r * pix * 2:.1f}, {width:.1f}, {pf}\n')
        print(f'{source_name}: {r * pix * 2:.1f}", {width:.1f} kpc,'
              f'{thresh_ans}')

        with open(results_csv, 'a') as f:
            f.write(result)
    return results_csv


def main():
    """Plot postage stamp images of LDR2 BL Lacs.
    """
    loop_through_sources()


if __name__ == '__main__':
    main()
