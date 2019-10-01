#!/usr/bin/env python3

"""Plot postage stamp images of LDR2 BL Lacs."""

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


def smallest_circle(sigma=4):
    """Measure the extent of blazars in LDR1.

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

    results_csv = f'{my_directory}/results/ldr2.bllacs.fsrqs.results.csv'
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

    big_sources = ['5BZBJ1202+4444', '5BZBJ1419+5423']

    for source_name, ra, dec, mosaic, rms, z in zip(
            df['Source name'], df['RA (J2000.0)'], df['Dec (J2000.0)'],
            df['Mosaic_ID'], df['Isl_rms'], df['Redshift']):

        field = f'{my_directory}/mosaics/{mosaic}-mosaic.fits'
        threshold = sigma * rms / 1000   # jansky
        save = f'{my_directory}/images/{source_name}.png'

        hdu = fits.open(field)[0]
        wcs = WCS(hdu.header, naxis=2)
        sky_position = SkyCoord(ra, dec, unit='deg')
        size = [3, 3] if source_name in big_sources else [2, 2]
        cutout = Cutout2D(np.squeeze(hdu.data), sky_position,
                          size=size * u.arcmin, wcs=wcs)

        d = cutout.data
        # https://docs.scipy.org/doc/numpy/reference/generated/numpy.copy.html
        copy = np.copy(d)
        d[d < threshold] = np.nan
        rows, cols = d.shape

        set_to_nan = []
        for r in range(rows):  # set to nan if surrounded by non nans
            for c in range(cols):
                try:
                    if (not np.isnan(d[r - 1, c - 1]) and
                            not np.isnan(d[r - 1, c]) and
                            not np.isnan(d[r - 1, c + 1]) and
                            not np.isnan(d[r, c - 1]) and
                            not np.isnan(d[r, c + 1]) and
                            not np.isnan(d[r + 1, c - 1]) and
                            not np.isnan(d[r + 1, c]) and
                            not np.isnan(d[r + 1, c + 1])):
                        set_to_nan.append((r, c))
                except IndexError:
                    continue

        for r, c in set_to_nan:
            d[r, c] = np.nan  # needs separate loop to avoid checkered pattern

        good_cells = []
        for r in range(rows):
            for c in range(cols):
                if not np.isnan(d[r, c]):
                    good_cells.append([r, c])

        # find distance between good_cell and all other good_cells
        max_distances, max_x1, max_x2, max_y1, max_y2 = [], [], [], [], []
        r = int(np.round(d.shape[0] / 2, 0))
        c = int(np.round(d.shape[1] / 2, 0))

        for (x1, y1) in good_cells:
            for (x2, y2) in good_cells:
                max_distances.append(np.sqrt((x1 - x2) ** 2 +
                                     (y1 - y2) ** 2))
                max_x1.append(x1)
                max_x2.append(x2)
                max_y1.append(y1)
                max_y2.append(y2)

        # not so simple - need the shortest of the lines that passes through
        # the centre
        d[r - 1:r + 1, c - 1:c + 1] = 0  # set centre to zero so we can see it

        max_distances = np.array(max_distances)
        my_max = np.max(max_distances)
        max_x1 = max_x1[max_distances.argmax()]
        max_x2 = max_x2[max_distances.argmax()]
        max_y1 = max_y1[max_distances.argmax()]
        max_y2 = max_y2[max_distances.argmax()]
        midpoint = ((max_x1 + max_x2) / 2, (max_y1 + max_y2) / 2)
        asec_max = my_max * 1.5  # 1.5" per pixel

        ax = plt.subplot(projection=wcs)
        plt.xlabel('Right ascension', fontsize=20, color='black')
        plt.ylabel('Declination', fontsize=20, color='black')
        ax.tick_params(axis='both', which='major', labelsize=20)
        plt.imshow(copy, vmin=0, vmax=np.nanmax(copy), origin='lower',
                   norm=DS9Normalize(stretch='arcsinh'), cmap='plasma_r')
        #           interpolation='gaussian',
        plt.plot([max_y1, max_y2], [max_x1, max_x2], color='black', alpha=1,
                 lw=2)
        beam = Circle((10, 10), radius=2, linestyle='dashed', color='grey',
                      lw=2)
        diffuse = Circle(midpoint, radius=my_max / 2, color='k', lw=2)
        ax.add_patch(beam)
        ax.add_patch(diffuse)
        cbar = plt.colorbar()
        cbar.set_label(r'Jy beam$^{-1}$', size=20)
        cbar.ax.tick_params(labelsize=20)
        plt.minorticks_on()
        plt.tick_params(which='minor', length=0)
        plt.contour(copy, levels=[threshold], origin='lower', colors='w')
        plt.savefig(save)
        plt.clf()

        dl, kpc = get_dl_and_kpc_per_asec(z=z)
        width = asec_max * kpc

        result = (f'{source_name},{ra},{dec},{rms * 1e3},{z}'
                  f',{asec_max:.1f},{width:.1f}\n')
        print(f'{source_name}: {asec_max:.1f}", {width:.2f} kpc')

        with open(results_csv, 'a') as f:
            f.write(result)
    return results_csv


def main():
    """Plot postage stamp images of LDR2 BL Lacs."""

    smallest_circle()


if __name__ == '__main__':
    main()
