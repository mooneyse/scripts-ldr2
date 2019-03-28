#!/usr/bin/env python3

'''Fit a Gaussian point spread function to a point source and subtract it from a source with diffuse emission.'''

import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib_scalebar.scalebar import ScaleBar
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy.wcs import WCS
from photutils.isophote import build_ellipse_model, Ellipse, EllipseGeometry
from scipy.ndimage.interpolation import map_coordinates, shift
import scipy.optimize as opt
from ds9norm import DS9Normalize

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '20 March 2019'

def get_df(filename, format, index):
    '''Create the data frame.'''
    if format is 'csv':
        df = pd.read_csv(filename)
    else:
        data = Table.read(filename, format=format)
        df = data.to_pandas()
    df.set_index(index, inplace=True)
    return df


def get_position(df, field='Bootes'):
    '''Look up the position of the blazar.'''
    df = df[df['Field'] == field]
    blazar_names = df.index.tolist()
    blazar_positions = []
    for ra, dec in zip(df['RA'], df['DEC']):
        blazar_positions.append([ra, dec])
    return blazar_names, blazar_positions


def nearest_point_source(df, position, s_code='S', flux_threshold=0.01,
                         distance_threshold=0.6, elongation_threshold=1.2):
    '''Find the nearest bright point source to a given blazar.'''
    pd.options.mode.chained_assignment = None  # disable warning, see https://stackoverflow.com/a/20627316/6386612
    df.reset_index(inplace=True)  # remove Source_id as index column
    df['S_Code'] = df['S_Code'].str.decode('utf-8')  # byte string encoded
    df_point_sources = df[(df['S_Code'] == s_code) & (df['Total_flux'] > flux_threshold)]
    blazar = SkyCoord(position[0], position[1], unit='deg')
    separations, elongations = [], []
    for ra, dec, major, minor in zip(df_point_sources['RA'], df_point_sources['DEC'], df_point_sources['Maj'], df_point_sources['Min']):
        point_source = SkyCoord(ra, dec, unit='deg')
        separation = blazar.separation(point_source)
        separations.append(separation.deg)
        elongations.append(major / minor)
    df_point_sources['Separation'] = separations
    df_point_sources['Elongation'] = elongations
    df_point_sources = df_point_sources[df_point_sources['Elongation'] < elongation_threshold]
    nearest = df_point_sources.loc[df_point_sources['Separation'].idxmin()]
    point_source_position = [nearest['RA'], nearest['DEC']]
    results = {'i': 'ILTJ' + str(nearest['Source_id']), 's': nearest['Separation'], 'f': nearest['Total_flux'] * 1000}
    print('{i} is {s:.2f} degrees away and has a total flux density of {f:.2f} mJy.'.format(**results))
    if results['s'] > distance_threshold:
        print('The point source is too far away.')
        sys.exit()
    return results['i'], point_source_position


def get_fits(filename):
    '''Open FITS image.'''
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header, naxis=2)
    return hdu, wcs


def get_data(position, hdu, wcs, size=[1, 1] * u.arcmin):
    '''Cut out the source from the FITS data.'''
    sky_position = SkyCoord(position[0], position[1], unit='deg')
    cutout = Cutout2D(np.squeeze(hdu.data), sky_position, size=size, wcs=wcs)
    data = cutout.data
    return data


def housekeeping(name, data):
    '''Some ad hoc adjustments to a few sources, made after visual inspection.'''
    if name == '5BZBJ1426+3404':
        print('Doing a little housekeeping on {}.'.format(name))
        data[:10, :10] = 0
    elif name == '5BZQJ1429+3529':
        print('Doing a little housekeeping on {}.'.format(name))
        data[:10, 20:] = 0
    elif name == '5BZQJ1435+3353':
        print('Doing a little housekeeping on {}.'.format(name))
        data[8:13, 28:33] = 0
    elif name == '5BZQJ1437+3618':
        print('Doing a little housekeeping on {}.'.format(name))
        data[1:10, :3] = 0
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
    return new_data  # this function was tested and worked as expected


def gaussian(xy, amplitude, x0, y0, sigma_x, sigma_y, theta, offset):
    '''Fit a two-dimensional Gaussian to the data.'''
    x, y = xy
    a = ((np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2))
    b = (-(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2))
    c = ((np.sin(theta) ** 2)/(2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2))
    g = (offset + amplitude * np.exp( -(a * ((x - x0) ** 2) + 2 * b * (x - x0) * (y - y0) + c * ((y - y0) **2 ))))
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
    '''Shift the data so that the index of the maximum value in the data matches up with that of the model. This is needed to ensure accurate subtraction.'''
    data_max = data.argmax()
    model_max = model.argmax()
    data_shape = data.shape
    model_shape = model.shape
    data_peak_row, data_peak_column = np.unravel_index(data_max, data_shape)
    model_peak_row, model_peak_column = np.unravel_index(model_max, model_shape)
    shifting = (model_peak_row - data_peak_row, model_peak_column - data_peak_column)
    shifted_data = shift(data, shift=shifting, cval=cval)
    numbers = {'d': (data_peak_row, data_peak_column), 's': shifting, 'm': (model_peak_row, model_peak_column)}
    print('The data were shifted from {d} by {s} to {m}.'.format(**numbers))
    return shifted_data


def get_noise_catalogue(df, blazar_name, new_data='', sigma=5):
    '''Calculate the noise in the blazar image.'''
    rms = df.loc[blazar_name, 'Isl_rms']
    five_sigma = [rms * sigma]
    return five_sigma


def jy_per_beam_to_jy(flux, beam_radius=6, sigma=1):
    '''Convert Jansky per beam measurements to Jansky measurements.'''
    bmaj = beam_radius * u.arcsec
    bmin = beam_radius * u.arcsec
    fwhm_to_sigma = 1 / np.sqrt(8 * np.log(2))
    beam_area = 2 * np.pi * bmaj * bmin * (fwhm_to_sigma ** 2)
    beam_area = beam_area.to(u.sr)
    flux_jy_per_beam = flux * u.Jy / u.sr
    flux_jy = flux_jy_per_beam * beam_area
    # BUG this is way off but the code looks right
    # see http://docs.astropy.org/en/stable/api/astropy.units.equivalencies.brightness_temperature.html
    return flux_jy


def diffuse_fraction(df, name, blazar, diffuse, threshold):
    '''Calculate what the fraction of emission is diffuse.'''
    blazar_catalogue = df.loc[name, 'Total_flux']
    blazar_image = blazar[blazar > threshold].sum()
    diffuse_image = diffuse[diffuse > threshold].sum()
    fraction = diffuse_image / blazar_image
    blazar_diffuse_catalogue = blazar_catalogue * fraction
    blazar_core_catalogue = blazar_catalogue * (1 - fraction)
    print('{0:.2%} of the emission is not from the core.'.format(fraction))
    results = {'total': blazar_catalogue * 1000, 'core': blazar_core_catalogue * 1000, 'diffuse': blazar_diffuse_catalogue * 1000}  # mJy
    print('The total flux density is {total:.2f} mJy, with the core contributing {core:.2f} mJy and the diffuse emission contributing {diffuse:.2f} mJy.'.format(**results))
    print('=', name+',', fraction)


def make_plot(position, data, title, rows=2, columns=7, origin='lower',
              vmin=0, vmax=1, zmin=0, zmax=1, axis='off', elev=24, azim=29,
              linewidth=0.25, pane_colour='white', line_colour='black',
              cmap='magma_r', projection='3d', stretch='arcsinh', levels='',
              plot='', layer='', contour_colours=['black', 'magenta'], pad=0.05,
              size='5%', orientation='horizontal', location='bottom',
              linestyles='dashed'):
    '''Plot the image on a grid.'''
    ax0 = plt.subplot(rows, columns, position)
    ax0.set_title(title)
    ax0.get_xaxis().set_ticks([])
    ax0.get_yaxis().set_ticks([])
    image_size = 60  # image is 1 * 1 arcminutes
    len_x, len_y = data.shape
    x0, x1 = len_x * 0.65, len_x * 0.9
    y = len_y * 0.1
    plt.text((x0 + x1) / 2, y * (1.4), str(int(image_size * 0.25)) + '"', horizontalalignment='center', verticalalignment='center', fontsize=12)
    plt.plot([x0, x1], [y, y], 'k-', lw=1.5)
    image = ax0.imshow(data, origin=origin, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=DS9Normalize(stretch=stretch))
    divider = make_axes_locatable(ax0)
    cax0 = divider.append_axes(location, size=size, pad=pad)
    plt.colorbar(image, cax=cax0, orientation=orientation)
    if plot is 'blazar':
        ax0.contour(data, levels=levels, origin=origin, colors=contour_colours[0])
    elif plot is 'diffuse':
        ax0.contour(layer, levels=levels, origin=origin, colors=contour_colours[0])
        ax0.contour(data, levels=levels, origin=origin, linestyles=linestyles, colors=contour_colours[1])
    len_x, len_y = data.shape
    range_x = range(len_x)
    range_y = range(len_y)
    ax1 = plt.subplot(rows, columns, int(position + columns), projection=projection)
    x, y = np.meshgrid(range_x, range_y)
    ax1.get_xaxis().set_ticks([])
    ax1.get_yaxis().set_ticks([])
    ax1.zaxis.set_tick_params(pad=10)  # ax1.set_zticks([]) to remove z labels
    ax1.xaxis.pane.fill = False
    ax1.yaxis.pane.fill = False
    ax1.zaxis.pane.fill = False
    ax1.xaxis.pane.set_edgecolor(pane_colour)
    ax1.yaxis.pane.set_edgecolor(pane_colour)
    ax1.zaxis.pane.set_edgecolor(pane_colour)
    ax1.grid(False)
    ax1.set_zlim(zmin, vmax)
    ax1.view_init(elev=elev, azim=azim)
    ax1.plot_surface(x, y, data, cmap=cmap, vmin=vmin, vmax=vmax, linewidth=linewidth, edgecolors=line_colour)


def main():
    '''Fit a Gaussian point spread function to a point source and subtract it from a source with diffuse emission.'''
    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=formatter_class)
    parser.add_argument('-i', '--image', required=False, type=str, default='/mnt/closet/ldr2-blazars/deep-fields/bootes-image.fits', help='FITS image of the field')
    parser.add_argument('-c', '--catalogue', required=False, type=str, default='/mnt/closet/ldr2-blazars/deep-fields/bootes.fits', help='FITS catalogue of the field')
    parser.add_argument('-d', '--csv', required=False, type=str, default='/mnt/closet/ldr2-blazars/deep-fields/bootes-lockman-hole-blazars.csv', help='CSV catalogue of the blazars')
    parser.add_argument('-o', '--output', required=False, type=str, default='/mnt/closet/ldr2-blazars/deep-fields/images/gaussian', help='Directory to save the plots')

    args = parser.parse_args()
    image = args.image
    catalogue = args.catalogue
    csv = args.csv
    output = args.output

    font = 'STIXGeneral'
    math_font = 'cm'
    font_size = 12
    figsize = (40, 20)
    bbox_inches = 'tight'
    testing = False
    new_size = 10

    df_blazars = get_df(csv, format='csv', index='Source name')
    blazar_names, blazar_positions = get_position(df_blazars)
    for i, (blazar_name, blazar_position) in enumerate(zip(blazar_names, blazar_positions)):
        if testing:
            if i != 0:  # do one at a time
                sys.exit()
        print('Analysing {} (blazar {} of {}).'.format(blazar_name, i + 1, len(blazar_names)))
        df_bootes = get_df(catalogue, format='fits', index='Source_id')
        point_source_id, point_source_position = nearest_point_source(df_bootes, blazar_position)
        hdu, wcs = get_fits(filename=image)
        blazar_data = get_data(position=blazar_position, hdu=hdu, wcs=wcs)
        blazar_data = housekeeping(blazar_name, blazar_data)
        point_source_data = get_data(position=point_source_position, hdu=hdu, wcs=wcs)
        blazar_regrid = regrid(blazar_data, new_size=new_size, normalise=False)  # peak and total values change with regridding
        point_source_regrid = regrid(point_source_data, new_size=new_size, normalise=False)
        model = make_model(point_source_regrid, sigma_x=4, sigma_y=4)
        point_source_residual = point_source_regrid - model
        scaled_model = (model * np.max(blazar_regrid) / np.max(point_source_regrid))
        blazar_shifted = match_peaks(blazar_regrid, scaled_model)
        blazar_residual = blazar_shifted - scaled_model
        blazar_regrid_back = regrid(blazar_shifted, new_size=1 / new_size, normalise=False)  # regrid the blazar and blazar residual data back to the native resolution
        blazar_residual_regrid_back = regrid(blazar_residual, new_size=1 / new_size, normalise=False)
        five_sigma = get_noise_catalogue(df_blazars, blazar_name)
        diffuse_fraction(df=df_blazars, name=blazar_name, blazar=blazar_regrid_back, diffuse=blazar_residual_regrid_back, threshold=five_sigma)
        savefig = output + '/' + blazar_name + '.png'
        matplotlib.rcParams['font.family'] = font
        matplotlib.rcParams['mathtext.fontset'] = math_font
        matplotlib.rcParams['font.size'] = font_size
        plt.figure(figsize=figsize)
        make_plot(position=1, data=point_source_regrid, title=point_source_id, vmax=np.max(point_source_regrid))
        make_plot(position=2, data=model, title=point_source_id + ' model', vmax=np.max(point_source_regrid))
        make_plot(position=3, data=point_source_residual, title=point_source_id + ' residual', vmax=np.max(point_source_regrid))
        make_plot(position=4, data=blazar_shifted, title=blazar_name, vmax=np.max(blazar_shifted))
        make_plot(position=5, data=blazar_residual, title=blazar_name + ' diffuse', vmax=np.max(blazar_shifted))
        make_plot(position=6, data=blazar_regrid_back, title=blazar_name, levels=five_sigma, plot='blazar', vmax=np.max(blazar_regrid_back))
        make_plot(position=7, data=blazar_residual_regrid_back, title=blazar_name + ' diffuse', levels=five_sigma, plot='diffuse', layer=blazar_regrid_back, vmax=np.max(blazar_regrid_back))
        if testing:
            plt.show()
        else:
            plt.savefig(savefig, bbox_inches=bbox_inches)


if __name__ == '__main__':
    main()
