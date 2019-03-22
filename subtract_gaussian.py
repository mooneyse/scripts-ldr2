#!/usr/bin/env python3

'''Fit a Gaussian point spread function to a point source and subtract it from
a source with diffuse emission.'''

import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.modeling.models import Gaussian2D
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.utils.data import get_pkg_data_filename
from photutils.datasets import make_noise_image
from photutils.isophote import build_ellipse_model, Ellipse, EllipseGeometry
from photutils import EllipticalAperture

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '20 March 2019'

def get_fits(filename):
    '''Open FITS file.'''
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header, naxis=2)
    return hdu, wcs


def get_data(position, hdu, wcs, size=[1, 1] * u.arcmin):
    '''Cut out the source from the FITS data.'''
    sky_position = SkyCoord(position[0], position[1], unit='deg')
    cutout = Cutout2D(np.squeeze(hdu.data), sky_position, size=size, wcs=wcs)
    data_normalised = cutout.data / np.max(cutout.data)
    return data_normalised


def make_model(data, x0=20, y0=19, sma=3, eps=0, pa=0):
    '''Fit a two-dimensional Gaussian to the data.'''
    geometry = EllipseGeometry(x0=x0, y0=y0, sma=sma, eps=eps, pa=pa)
    ellipse = Ellipse(data, geometry)
    iso_list = ellipse.fit_image()
    model = build_ellipse_model(data.shape, iso_list)
    return model


def make_plot(position, data, title, rows=2, columns=3, origin='lower', vmin=0,
              vmax=1, axis='off'):
    '''Plot the image on a grid.'''
    ax = plt.subplot(int(str(rows) + str(columns) + str(position)))
    ax.imshow(data, origin=origin, vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.axis(axis)


def main():
    '''Fit a Gaussian point spread function to a point source and subtract it
    from a source with diffuse emission.'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-f',
                        '--filename',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/deep-fields/bootes-image.fits',
                        help='FITS file of the field')

    parser.add_argument('-p',
                        '--point',
                        required=False,
                        nargs='+',
                        default=[216.71718, 33.97577],
                        help='FITS file containing a point source')

    parser.add_argument('-b',
                        '--blazar',
                        required=False,
                        nargs='+',
                        default=[216.5321226, 34.07423184],
                        help='FITS file containing a diffuse source')

    args = parser.parse_args()
    filename = args.filename
    point_source_position = args.point
    blazar_position = args.blazar

    hdu, wcs = get_fits(filename=filename)
    blazar_data = get_data(position=blazar_position, hdu=hdu, wcs=wcs)
    point_source_data = get_data(position=point_source_position, hdu=hdu, wcs=wcs)
    model = make_model(data=point_source_data)
    point_source_residual = point_source_data - model
    blazar_residual = blazar_data - model

    make_plot(position=1, data=point_source_data, title='Point source')
    make_plot(position=2, data=model, title='Point source model')
    make_plot(position=3, data=point_source_residual, title='Point source residual')
    make_plot(position=4, data=blazar_data, title='5BZBJ1426+3404')
    make_plot(position=6, data=blazar_residual, title='5BZBJ1426+3404 residual')

    plt.show()


if __name__ == '__main__':
    main()
