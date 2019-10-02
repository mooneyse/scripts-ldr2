#!/usr/bin/env python3

"""Blazar LDR2 analysis.
"""

import os
import argparse
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib as mpl
# mpl.use('Agg')
from astropy.wcs import WCS
from astropy.io import fits
from ds9norm import DS9Normalize

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '18 August 2019'


def kpc_per_asec(z, H0=70, WM=0.26, WV=0.74, n=1000):
    """Ned Wright's cosmology calculator.

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
    h = H0 / 100
    WR = 4.165E-5 / (h * h)
    WK = 1 - WM - WR - WV
    az = 1 / (1 + 1 * z)
    age, DCMR = 0, 0

    for i in range(n):
        a = az * (i + 0.5) / n
        adot = math.sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
        age = age + 1 / adot

    for i in range(n):
        a = az + (1 - az) * (i + 0.5) / n
        adot = math.sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
        DCMR = DCMR + 1 / (a * adot)

    DCMR = (1 - az) * DCMR / n
    x = math.sqrt(abs(WK)) * DCMR

    if x > 0.1:
        if WK > 0:
            ratio = 0.5 * (math.exp(x) - math.exp(-x)) / x
        else:
            ratio = math.sin(x) / x
    else:
        y = x * x
        if WK < 0:
            y = -y
        ratio = 1 + y / 6 + y * y / 120

    return (299792.458 / H0) * az * ratio * (DCMR / 206.264806)


def main():
    """Blazar LDR2 analysis.

    Parameters
    ----------

    Returns
    -------
    None
    """

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)
    parser.add_argument('-c',
                        '--csv',
                        required=False,
                        type=str,
                        default=('/mnt/closet/ldr2-blazars/catalogues/' +
                                 'bl_lacs_with_z.csv'),
                        help='CSV catalogue of the blazars with redshifts')

    args = parser.parse_args()
    csv = args.csv
    directory = '/mnt/closet/ldr2-blazars/images/fits/bzb'
    save = '/mnt/closet/ldr2-blazars/images/png/bzb-for-fra/with-z-with-circle'
    zoom = 8  # factor to zoom in by
    sigma = 4  # signifiance of the contour

    # already split sample into those with and without a reliable redshift;
    # focusing on the sample with redshifts

    # for each source in the sample get the redshift
    # use the redshift to get the kiloparsecs per pixel
    df = pd.read_csv(csv)
    df = df[df['Name'] == '5BZB J0942+2844']  # filter one source

    # draw a 30 kpc diameter circle on the radio images at the peak flux
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
    mpl.rcParams['axes.edgecolor'] = 'black'

    for name, ra, dec, pointing, z in zip(df['Name'], df['RA'], df['Dec'],
                                          df['Pointing'], df['Redshift']):
        if isinstance(name, str):  # nans are float this so remove them
            name = name.replace(' ', '')
            my_path = f'{directory}/{name}-{pointing}-DR2-mJy.fits'

            if os.path.isfile(path=my_path):  # see if we have the fits cutout
                kpc_sec = kpc_per_asec(z=z)  # with 1.5 pixels per arcsecond
                kpc = (1 / kpc_sec) * 30
                kpc = kpc / 60 / 60  # degrees

                # get the data
                hdu = fits.open(my_path)[0]
                data = np.squeeze(hdu.data)

                # zoom in as specified
                x, _ = data.shape
                a, b = int(x / 2 * (1 - 1 / zoom)), int(x / 2 * (1 + 1 / zoom))
                data = data[a:b, a:b]
                wcs = WCS(hdu.header, naxis=2)[a:b, a:b]

                # get the significance contour
                # using the residuals overestimates the noise
                residual = f'/mnt/closet/ldr2-blazars/images/noise/{name}.fits'
                try:
                    rms = np.std(np.squeeze(fits.open(residual)[0].data))
                except FileNotFoundError:  # do not have the residual
                    rms = 0.6  # an estimate
                level = rms * sigma

                # make the plot
                ax = plt.subplot(projection=wcs)
                ax.tick_params(direction='in', length=6, width=2)
                print(kpc*60*60)
                galaxy = Circle((ra, dec), kpc/2, edgecolor='k',
                                facecolor='none',
                                transform=ax.get_transform('icrs'))
                # add the 6" (4 pixel) lofar beam to the images
                beam = Circle((4, 4), 2, edgecolor='yellow', facecolor='none',
                              linestyle='--')
                ax.add_patch(galaxy)
                # ax.set_autoscale_on(False)
                ax.add_patch(beam)
                plt.xlabel('Right ascension', fontsize=20, color='black')
                plt.ylabel('Declination', fontsize=20, color='black')
                ax.tick_params(axis='both', which='major', labelsize=20)
                plt.imshow(data, origin='lower', vmin=0, vmax=np.max(data),
                           norm=DS9Normalize(stretch='arcsinh'))
                ax.contour(data, levels=[level], colors=['magenta'])

                # panstarrs contour
                panstarrs = '/home/sean/Downloads/cutout_rings.v3.skycell.1962.015.stk.i.unconv.fits'
                p_hdu = fits.open(panstarrs)[0]
                p_data = p_hdu.data
                p_wcs = WCS(p_hdu.header)
                ax.contour(p_data, levels=[845], colors='w', transform=ax.get_transform(p_wcs))  # 845 is 5 sigma from a manual calculation
                plt.xlim(0, data.shape[1] - 1)
                plt.ylim(0, data.shape[0] - 1)
                plt.title(f'{name}, z = {z}, 1" = {kpc_sec:.1f} kpc',
                          fontsize=20)
                plt.savefig(f'{save}/{name}.png')
                plt.clf()
                # break

    # https://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=4.40375%2C14.8505277777778&filter=g&filter=r&filter=i&filter=z&filter=y&filetypes=stack&size=240&output_size=0&verbose=0&autoscale=99.500000&catlist=

    # use panstarrs to mark the position of optical sources in the field

    # Elliptical galaxies have very little gas and dust. Since stars form from gas, little star formation occurs in elliptical galaxies. Most of their stars are old and red.

    # inspect images to check for source confusion, splitting the sample into
    # bl lacs with and without source confusion

    # for those without source confusion, add the smallest circle possible to
    # encompass the 5 sigma emission

    # report the size of the circle to a csv

    # split the bl lac sample into those with and without emission past 30 kpc

    # remove the tev blazars for the moment

    # ask Tim what the ldr2 catalogue and publication schedule is

    # organise the next meeting for the end of august and include the work I
    # have completed

    # check to see if the sample fra sent on of elliptical galaxies with z <
    # 0.3 and no nvss counterparts have lofar detections

    # update the wiki, including two images: one with and one without two
    # sources in the fov (i.e. give an example of what we are doing)

    # text in the wiki will become the paper so add details of the method here
    # as it progresses.

    # TODO instead of a 30 kpc circle, download the panstarrs images for each
    # source
    # plot radio contours over them
    # fit the smallest circle to the optical contours
    # report the size of the optical smallest circle to the csv
    # see how many sources have radio emission beyond the optical emission

    # TODO introduce a flag into the csv to indicate whether an image contains
    # artefacts (look at the residual images to see if there are outliers in
    # it)

    # TODO plot first contours on all images (especially for 1617+4106)


if __name__ == '__main__':
    main()
