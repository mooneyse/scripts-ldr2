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
from scipy.ndimage.interpolation import map_coordinates, shift
import scipy.optimize as opt
from ds9norm import DS9Normalize

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
    data = cutout.data
    return data


def regrid(data, new_size=10, normalise=True):
    '''Map the data onto a larger array.'''
    l, w = data.shape
    new_l = l * new_size
    new_w = w * new_size
    new_dimensions = []
    for old_shape, new_shape in zip(data.shape, (new_l, new_w)):
        new_dimensions.append(np.linspace(0, old_shape - 1, new_shape))

    coordinates = np.meshgrid(*new_dimensions, indexing='ij')
    new_data = map_coordinates(data, coordinates)
    new_data = new_data / np.max(new_data) if normalise else new_data
    return new_data


def gaussian(xy, amplitude, x0, y0, sigma_x, sigma_y, theta, offset):
    '''Fit a two-dimensional Gaussian to the data.'''
    x, y = xy
    a = ((np.cos(theta) ** 2) / (2 * sigma_x ** 2) +
    (np.sin(theta) ** 2) / (2 * sigma_y ** 2))
    b = (-(np.sin(2 * theta)) / (4 * sigma_x ** 2) +
    (np.sin(2 * theta)) / (4 * sigma_y ** 2))
    c = ((np.sin(theta) ** 2)/(2 * sigma_x ** 2) +
    (np.cos(theta) ** 2) / (2 * sigma_y ** 2))
    g = (offset + amplitude * np.exp( -(a * ((x - x0) ** 2) +
    2 * b * (x - x0) * (y - y0) + c * ((y - y0) **2 ))))
    g_ravel = g.ravel()
    return g_ravel


def make_model(data, amplitude=1, sigma_x=30, sigma_y=30, theta=0, offset=0):
    '''Fit a model to the data.'''
    len_x, len_y = data.shape
    range_x = range(len_x)
    range_y = range(len_y)
    x, y = np.meshgrid(range_x, range_y)
    x0 = len_x / 2
    y0 = len_y / 2
    data_ravel = data.ravel()
    p0 = (amplitude, x0, y0, sigma_x, sigma_y, theta, offset)
    popt, pcov = opt.curve_fit(gaussian, (x, y), data_ravel, p0=p0)
    model = gaussian((x, y), *popt)
    model_reshape = model.reshape(len_x, len_y)
    return model_reshape


def match_peaks(data, model, cval=0):
    '''Shift the data so that the index of the maximum value in the data
    matches up with that of the model. This is needed to ensure accurate
    subtraction.'''
    data_max = data.argmax()
    model_max = model.argmax()
    data_shape = data.shape
    model_shape = model.shape
    data_peak_row, data_peak_column = np.unravel_index(data_max, data_shape)
    model_peak_row, model_peak_column = np.unravel_index(model_max, model_shape)
    shifting = (model_peak_row - data_peak_row, model_peak_column - data_peak_column)
    shifted_data = shift(data, shift=shifting, cval=cval)
    numbers = {'d': (data_peak_row, data_peak_column), 's': shifting,
               'm': (model_peak_row, model_peak_column)}
    print('The data were shifted from {d} by {s} to match {m}.'.format(**numbers))
    return shifted_data


def get_noise(original_data, new_data, x0=0, x1=10, y0=0, y1=10, sigma=5):
    '''Calculate the noise in the blazar image.'''
    maximum = np.max(original_data)
    deviation = np.std(original_data[x0:x1, y0:y1])
    snr = maximum / deviation
    print('The signal to noise ratio is {:.3f}'.format(snr))
    threshold = [(np.max(new_data) / snr) * sigma]  # list for contouring
    return threshold


def diffuse_fraction(data, residual, threshold):
    '''Calculate what the fraction of emission is diffuse.'''
    blazar_total = data[data > threshold].sum()
    diffuse_total = residual[residual > threshold].sum()
    fraction = diffuse_total / blazar_total
    print('{0:.2%} of the blazar emission is not from the core.'.format(fraction))
    return fraction


def make_plot(position, data, title, rows=2, columns=5, origin='lower',
              vmin=0, vmax=1, zmin=0, zmax=1, axis='off', elev=28, azim=-61,
              linewidth=0.25, pane_colour='white', line_colour='black',
              cmap='magma_r', projection='3d', stretch='arcsinh', alpha=1,
              levels='', source=''):
    '''Plot the image on a grid.'''
    ax0 = plt.subplot(rows, columns, position)
    ax0.set_title(title)
    ax0.get_xaxis().set_ticks([])
    ax0.get_yaxis().set_ticks([])
    ax0.imshow(data, origin=origin, vmin=vmin, vmax=vmax, cmap=cmap,
               norm=DS9Normalize(stretch=stretch))
    if source is not '':
        ax0.contour(source, levels=levels, origin=origin, alpha=alpha, linestyles='dashed')
    if levels is not '':
        ax0.contour(data, levels=levels, origin=origin, alpha=alpha)

    len_x, len_y = data.shape
    range_x = range(len_x)
    range_y = range(len_y)

    ax1 = plt.subplot(rows, columns, int(position + columns), projection=projection)
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
    ax1.plot_surface(x, y, data, cmap=cmap, vmin=vmin, vmax=vmax,
                     linewidth=linewidth, edgecolors=line_colour)


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
    figsize= (16, 8)
    savefig = '/home/sean/Downloads/images/gaussian.png'

    # to get an automatic point source, search for the nearest source to a blazar
    # that was fit with a single gaussian and has a flux above a certain amount
    # the blazar name can also be read from the csv file, as can the island rms
    # i still need a way to convert these values to flux density

    hdu, wcs = get_fits(filename=filename)
    blazar_data = get_data(position=blazar_position, hdu=hdu, wcs=wcs)
    point_source_data = get_data(position=point_source_position, hdu=hdu, wcs=wcs)
    blazar_regrid = regrid(blazar_data)
    point_source_regrid =regrid(point_source_data)
    model = make_model(point_source_regrid)
    blazar_shifted = match_peaks(blazar_regrid, model)
    point_source_residual = point_source_regrid - model
    blazar_residual = blazar_shifted - model
    threshold = get_noise(blazar_data, blazar_shifted, x1=12, y1=12)  # box bounds determined by inspection
    diffuse_emission = diffuse_fraction(blazar_shifted, blazar_residual, threshold)

    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    matplotlib.rcParams['mathtext.fontset'] = 'cm'
    plt.figure(figsize=figsize)
    make_plot(position=1, data=point_source_regrid, title='Point source')
    make_plot(position=2, data=model, title='Point source model')
    make_plot(position=3, data=point_source_residual, title='Point source residual')
    make_plot(position=4, data=blazar_shifted, title='5BZBJ1426+3404', levels=threshold)
    make_plot(position=5, data=blazar_residual, title='5BZBJ1426+3404 residual', levels=threshold, source=blazar_shifted)
    # plt.show()
    plt.savefig(savefig, bbox_inches='tight')


if __name__ == '__main__':
    main()
