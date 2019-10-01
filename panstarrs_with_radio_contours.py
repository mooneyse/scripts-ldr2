#!/usr/bin/env python3

"""Plot Pan-STARRS data.
"""

import matplotlib as mpl
mpl.use('Agg')
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from ds9norm import DS9Normalize
from math import sqrt, exp, sin
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from skimage.measure import label
import smallestenclosingcircle

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '06 September 2019'


def nearest_to_centre(my_arr, percent):
    """Given a two dimensional array, return the value of the pixel nearest to
    the centre that is non-zero.
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
    else:
        R = int(round(my_arr.shape[0] / 2, 0))
        C = int(round(my_arr.shape[0] / 2, 0))

        dist = {}
        for r in range(my_arr.shape[0]):
            for c in range(my_arr.shape[1]):
                if my_arr[r, c] != 0:
                    dist[my_arr[r, c]] = np.sqrt((R - r) ** 2 + (C - c) ** 2)
        return [min(dist.items(), key=operator.itemgetter(1))[0]]


def get_dl_and_kpc_per_asec(z, H0=70, WM=0.26, WV=0.74):
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
        The luminosity distance in metres.
    float
        The kpc per arcsecond at the given redshift.
    """
    WR = 0.        # Omega(radiation)
    WK = 0.        # Omega curvaturve = 1-Omega(total)
    c = 299792.458  # velocity of light in km/sec
    DTT = 0.5      # time from z to now in units of 1/H0
    age = 0.5      # age of Universe in units of 1/H0
    zage = 0.1     # age of Universe at redshift z in units of 1/H0
    DCMR = 0.0     # comoving radial distance in units of c/H0
    DA = 0.0       # angular size distance
    DA_Mpc = 0.0
    kpc_DA = 0.0
    DL = 0.0       # luminosity distance
    DL_Mpc = 0.0
    a = 1.0        # 1/(1+z), the scale factor of the Universe
    az = 0.5       # 1/(1+z(object))
    h = H0/100.
    WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
    WK = 1-WM-WR-WV
    az = 1.0/(1+1.0*z)
    age = 0.
    n = 1000         # number of points in integrals
    for i in range(n):
        a = az*(i+0.5)/n
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        age = age + 1./adot

    zage = az*age/n
    DTT = 0.0
    DCMR = 0.0

    for i in range(n):
        a = az+(1-az)*(i+0.5)/n
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        DTT = DTT + 1./adot
        DCMR = DCMR + 1./(a*adot)

    DTT = (1.-az)*DTT/n
    DCMR = (1.-az)*DCMR/n
    age = DTT+zage

    ratio = 1.00
    x = sqrt(abs(WK))*DCMR
    if x > 0.1:
        if WK > 0:
            ratio = 0.5*(exp(x)-exp(-x))/x
        else:
            ratio = sin(x)/x
    else:
        y = x*x
        if WK < 0:
            y = -y
        ratio = 1. + y/6. + y*y/120.
    DCMT = ratio*DCMR
    DA = az*DCMT
    DA_Mpc = (c/H0)*DA
    kpc_DA = DA_Mpc/206.264806
    DL = DA/(az*az)
    DL_Mpc = (c/H0)*DL

    ratio = 1.00
    x = sqrt(abs(WK))*DCMR
    if x > 0.1:
        if WK > 0:
            ratio = (0.125*(exp(2.*x)-exp(-2.*x))-x/2.)/(x*x*x/3.)
        else:
            ratio = (x/2. - sin(2.*x)/4.)/(x*x*x/3.)
    else:
        y = x*x
        if WK < 0:
            y = -y
        ratio = 1. + y/5. + (2./105.)*y*y
    return DL_Mpc * 3.086e22, kpc_DA


def optical(sigma=4):
    """Plot Pan-STARRS data.

    Parameters
    ----------
    sigma : float or int
        The threshold of the significance to set the mask, as a factor of the
        local RMS. The default is 4.

    Returns
    -------
    string
        The name of the CSV containing the results.
    """
    my_directory = '/data5/sean/ldr2'
    results_csv = f'{my_directory}/results/ldr2.csv'
    df = pd.read_csv(f'{my_directory}/catalogues/final.csv')
    result_header = ('Name,RA,Dec,RMS (uJy),Redshift,Width ("),Width (kpc)\n')

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

    for source_name, ra, dec, mosaic, rms, z in zip(df['Source name'],
                                                    df['RA (J2000.0)'],
                                                    df['Dec (J2000.0)'],
                                                    df['Mosaic_ID'],
                                                    df['Isl_rms'],
                                                    df['Redshift']):

        sky_position = SkyCoord(ra, dec, unit='deg')
        if source_name == '5BZBJ1202+4444':
            size = [3, 3] * u.arcmin
        elif source_name == '5BZBJ1419+5423':
            size = [4, 4] * u.arcmin
        else:
            size = [2, 2] * u.arcmin

        hdu1 = fits.open(f'{my_directory}/panstarrs/{source_name}.i.fits')[0]
        wcs1 = WCS(hdu1.header)
        panstarrs = Cutout2D(np.squeeze(hdu1.data), sky_position, size=size,
                             wcs=wcs1).data

        hdu2 = fits.open(f'{my_directory}/mosaics/{mosaic}-mosaic.fits')[0]
        wcs2 = WCS(hdu2.header, naxis=2)
        ldr2 = Cutout2D(np.squeeze(hdu2.data), sky_position, size=size,
                        wcs=wcs2).data

        levels = [level * rms / 1000 for level in [4, 8, 16, 32]]
        save = f'{my_directory}/images/panstarrs-{source_name}.png'

        ax = plt.subplot(projection=wcs1)
        plt.xlabel('Right ascension', fontsize=20, color='black')
        plt.ylabel('Declination', fontsize=20, color='black')
        ax.tick_params(axis='both', which='major', labelsize=20)
        plt.imshow(panstarrs, vmin=0, vmax=8000,
                   origin='lower', norm=DS9Normalize(stretch='arcsinh'),
                   cmap='plasma_r')  # interpolation='gaussian'
        # plt.plot([panstarrs.shape[0] * 0.2, panstarrs.shape[0] * 0.2],
        #          [panstarrs.shape[1] * 0.25, panstarrs.shape[1] * 0.5],)
        cbar = plt.colorbar()
        cbar.set_label(r'Excess counts', size=20)
        cbar.ax.tick_params(labelsize=20)
        plt.minorticks_on()
        plt.tick_params(which='minor', length=0)
        plt.contour(ldr2, levels=levels, origin='lower',
                    colors=['red', 'yellow', 'blue', 'purple'])

        plt.savefig(save)
        plt.clf()

        dl, kpc = get_dl_and_kpc_per_asec(z=z)
        width = r * kpc

        result = (f'{source_name},{ra},{dec},{rms * 1e3},{z}'
                  f',{r:.1f},{width:.1f}\n')
        print(f'{source_name}: {r:.1f}", {width:.1f} kpc')

        with open(results_csv, 'a') as f:
            f.write(result)
    return results_csv


def main():
    """Plot Pan-STARRS data.
    """

    optical()


if __name__ == '__main__':
    main()
