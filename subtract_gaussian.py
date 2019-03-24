#!/usr/bin/env python3

'''Fit a Gaussian point spread function to a point source and subtract it from
a source with diffuse emission.'''

import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from photutils.isophote import build_ellipse_model, Ellipse, EllipseGeometry

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
    data_max = np.max(cutout.data)
    data_normalised = cutout.data / data_max
    return data_normalised, data_max


def make_model(data, x0=20, y0=19, sma=3, eps=0, pa=0):
    '''Fit a two-dimensional Gaussian to the data.'''
    geometry = EllipseGeometry(x0=x0, y0=y0, sma=sma, eps=eps, pa=pa)
    ellipse = Ellipse(data, geometry)
    iso_list = ellipse.fit_image()
    model = build_ellipse_model(data.shape, iso_list)
    return model


def make_plot(position, data, title, rows=2, columns=5, origin='lower',
              vmin=-1, vmax=1, zmin=-0.5, zmax=1, axis='off', elev=28, azim=-61,
              linewidth=0.25, pane_colour='white', line_colour='black',
              cmap='RdGy_r'):
    '''Plot the image on a grid.'''

    ax0 = plt.subplot(rows, columns, position)
    ax0.set_title(title)
    ax0.get_xaxis().set_ticks([])
    ax0.get_yaxis().set_ticks([])
    ax0.imshow(data, origin=origin, vmin=vmin, vmax=vmax, cmap=cmap)

    len_x, len_y = data.shape
    range_x = range(len_x)
    range_y = range(len_y)

    ax1 = plt.subplot(rows, columns, int(position + columns), projection='3d')
    x, y = np.meshgrid(range_x, range_y)
    ax1.get_xaxis().set_ticks([])
    ax1.get_yaxis().set_ticks([])
    ax1.set_zticks([])
    ax1.xaxis.pane.fill = False
    ax1.yaxis.pane.fill = False
    ax1.zaxis.pane.fill = False
    ax1.xaxis.pane.set_edgecolor(pane_colour)
    ax1.yaxis.pane.set_edgecolor(pane_colour)
    ax1.zaxis.pane.set_edgecolor(pane_colour)
    ax1.grid(False)
    ax1.set_zlim(zmin, zmax)
    ax1.view_init(elev=elev, azim=azim)
    ax1.plot_surface(x, y, data, cmap=cmap, vmin=vmin, vmax=vmax, alpha =0.9,
                     linewidth=linewidth, edgecolors=line_colour)
    ax1.contourf(x, y, data, zdir='z', offset=zmin, vmin=vmin, vmax=vmax,
                 cmap=cmap)


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
    figsize= (20, 10)
    savefig = '/home/sean/Downloads/images/gaussian.png'

    # to get an automatic point source, search for the nearest source to a blazar
    # that was fit with a single gaussian and has a flux above a certain amount
    # the blazar name can also be read from the csv file, as can the island rms
    # i still need a way to convert these values to flux density

    hdu, wcs = get_fits(filename=filename)
    blazar_data, blazar_max = get_data(position=blazar_position, hdu=hdu, wcs=wcs)
    point_source_data, point_source_max = get_data(position=point_source_position, hdu=hdu, wcs=wcs)
    model = make_model(data=point_source_data)
    point_source_residual = point_source_data - model
    blazar_residual = blazar_data - model

    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    matplotlib.rcParams['mathtext.fontset'] = 'cm'
    plt.figure(figsize=figsize)
    make_plot(position=1, data=point_source_data, title='Point source')
    make_plot(position=2, data=model, title='Point source model')
    make_plot(position=3, data=point_source_residual, title='Point source residual')
    make_plot(position=4, data=blazar_data, title='5BZBJ1426+3404')
    make_plot(position=5, data=blazar_residual, title='5BZBJ1426+3404 residual')
    plt.show()
    # plt.savefig(savefig, bbox_inches='tight')


    noise = 6.817118e-5
    sigma = 5
    threshold = noise * sigma

    # un-normalise the data
    blazar_unnormalised = blazar_residual * blazar_max
    # plt.imshow(blazar_unnormalised)
    # plt.show()
    # print(np.sum(blazar_unnormalised), 'thresh')
    # print(blazar_unnormalised[blazar_unnormalised > threshold].sum())



if __name__ == '__main__':
    main()
