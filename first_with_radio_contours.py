#!/usr/bin/env python3

"""Plot FIRST data.
"""

import matplotlib as mpl
# mpl.use('Agg')
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


def ghz(sigma=4):
    """Plot FIRST data.

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

    for source_name, ra, dec, mosaic, rms, z, f_rms in zip(df['Source name'],
                                                           df['RA (J2000.0)'],
                                                           df['Dec (J2000.0)'],
                                                           df['Mosaic_ID'],
                                                           df['Isl_rms'],
                                                           df['Redshift'],
                                                           df['FIRST RMS (mJy)']):

        sky_position = SkyCoord(ra, dec, unit='deg')
        if source_name == '5BZBJ1202+4444':
            size = [3, 3] * u.arcmin
        elif source_name == '5BZBJ1419+5423':
            size = [4, 4] * u.arcmin
        else:
            size = [2, 2] * u.arcmin

        first = f'{my_directory}/first/{source_name}.fits'
        f_hdu = fits.open(first)[0]
        f_wcs = WCS(f_hdu.header, naxis=2)
        f_cutout = Cutout2D(np.squeeze(f_hdu.data), sky_position, size=size,
                            wcs=f_wcs)

        ldr2 = f'{my_directory}/mosaics/{mosaic}-mosaic.fits'
        l_hdu = fits.open(ldr2)[0]
        l_wcs = WCS(l_hdu.header, naxis=2)
        l_cutout = Cutout2D(np.squeeze(l_hdu.data), sky_position, size=size,
                            wcs=l_wcs)

        levels = [level * rms / 1000 for level in [4, 8, 16, 32]]
        colors = ['#118ab2', '#06d6a0', '#ffd166', '#ef476f']

        f_level = f_rms * 4 / 1000
        fmt = {}
        fmt[f_level] = 'FIRST'

        ax = plt.subplot(projection=f_cutout.wcs)
        im = ax.imshow(f_cutout.data, vmin=0, vmax=np.max(f_cutout.data),
                       cmap='Greys', origin='lower',
                       norm=DS9Normalize(stretch='arcsinh'))
        cs = ax.contour(f_cutout.data, levels=[f_level], origin='lower',
                        colors=['k'])
        # ax.clabel(cs), fmt=fmt)
        # interpolation='gaussian'

        cbar = plt.colorbar(im)
        cbar.set_label(r'Jy beam$^{-1}$', size=20)
        cbar.ax.tick_params(labelsize=20)

        ax.contour(l_cutout.data, transform=ax.get_transform(l_cutout.wcs),
                   levels=levels, origin='lower', colors=colors)

        plt.xlabel('Right ascension', fontsize=20, color='black')
        plt.ylabel('Declination', fontsize=20, color='black')
        ax.tick_params(axis='both', which='major', labelsize=20)
        plt.minorticks_on()
        ax.tick_params(which='minor', length=0)

        plt.xlim(0, f_cutout.data.shape[0])
        plt.ylim(0, f_cutout.data.shape[1])

        plt.show()
        # return
        # save = f'{my_directory}/images/panstarrs-{source_name}.png'
        # plt.savefig(save)
        # plt.clf()

        #     dl, kpc = get_dl_and_kpc_per_asec(z=z)
        #     width = r * kpc


def main():
    """Plot FIRST data.
    """

    ghz()


if __name__ == '__main__':
    main()
