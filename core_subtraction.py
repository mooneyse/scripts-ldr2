#!/usr/bin/env python3

'''Fit a Gaussian point spread function to a point source and subtract it from
a source with diffuse emission.'''

import warnings
warnings.filterwarnings('ignore')  # supress warnings

import aplpy
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib_scalebar.scalebar import ScaleBar
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy.wcs import WCS
# from photutils.isophote import build_ellipse_model, Ellipse, EllipseGeometry
from scipy.ndimage.interpolation import map_coordinates, shift
import scipy.optimize as opt
from ds9norm import DS9Normalize

# def get_df(filename, format, index):
#     """Create the data frame.
#     """
#     if format == 'csv':
#         df = pd.read_csv(filename)
#     else:
#         data = Table.read(filename, format=format)
#         df = data.to_pandas()
#     df.set_index(index, inplace=True)
#     return df


# def get_position(df, cat_dir):
#     """Look up the position of the blazar.
#     """
#     blazar_names = df.index.tolist()
#     blazar_positions, catalogues, fits_images = [], [], []
#     for ra, dec, field in zip(df['BZCAT RA'], df['BZCAT Dec'],
#                               df['Mosaic_ID']):
#         blazar_positions.append([ra, dec])
#         field = field.lower().replace(' ', '.')
#         catalogues.append(f'{cat_dir}/ldr2-point-sources-near-bllacs.csv')
#         fits_images.append(f'/data5/sean/ldr2/mosaics/{field}-mosaic.fits')
#     return blazar_names, blazar_positions, catalogues, fits_images


# def new_nearest_point_source(csv, bllac):
#     df = pd.read_csv('/data5/sean/ldr2/catalogues/bright.near.points.csv')
#     df = df[df['BZCAT name'] == bllac]


def nearest_point_source(df, position, s_code='S', flux_threshold=0.01,
                         distance_threshold=1.5, elongation_threshold=1.2):
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
        return False, False
    return results['i'], point_source_position


def get_fits(filename):
    '''Open FITS image.'''
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header, naxis=2)
    return hdu, wcs


def get_data(position, hdu, wcs, size=[2, 2] * u.arcmin):
    '''Cut out the source from the FITS data.'''
    sky_position = SkyCoord(position[0], position[1], unit='deg')
    cutout = Cutout2D(np.squeeze(hdu.data), sky_position, size=size, wcs=wcs)
    data = cutout.data
    return data


def housekeeping(name, data):
    '''Some ad hoc adjustments to a few sources, made after visual inspection.
    '''
    if name == '5BZQJ1422+3223':
        print('Doing a little housekeeping on {}.'.format(name))
        data[42:47, 55:60] = 0
        data[6:10, 74:79] = 0

    if name == '5BZBJ1426+3404':
        print('Doing a little housekeeping on {}.'.format(name))
        data[20:28, 23:32] = 0  # data[:10, :10] = 0
        data[68:75, 33:42] = 0

    if name == '5BZQJ1429+3529':
        print('Doing a little housekeeping on {}.'.format(name))
        data[40:48, 10:17] = 0
        data[10:26, 40:60] = 0  # data[:10, 20:] = 0

    if name == '5BZQJ1435+3353':
        print('Doing a little housekeeping on {}.'.format(name))
        data[38:45, 5:11] = 0  # data[8:13, 28:33] = 0
        data[73:79, 54:60] = 0
        data[27:33, 48:53] = 0
        data[16:22, 64:71] = 0

    if name == '5BZQJ1437+3618':
        print('Doing a little housekeeping on {}.'.format(name))
        data[0:31, 4:26] = 0  # data[1:10, :3] = 0

    if name == '5BZQJ1437+3519':
        print('Doing a little housekeeping on {}.'.format(name))
        # data[:, 60:] = 0

    if name == '5BZBJ1558+5625':
        print('Doing a little housekeeping on {}.'.format(name))
        data[65:70, 36:43] = 0

    if name == '5BZBJ1605+5421':
        print('Doing a little housekeeping on {}.'.format(name))
        data[53:64, 16:27] = 0
        data[72:79, 35:45] = 0
        data[76:79, 76:80] = 0
        data[48:56, 75:80] = 0
        data[6:11, 70:76] = 0
        data[24:31, 42:48] = 0

    if name == '5BZQJ1606+5405':
        print('Doing a little housekeeping on {}.'.format(name))
        data[57:68, 69:79] = 0
        data[25:30, 11:17] = 0

    if name == '5BZQJ1608+5613':
        print('Doing a little housekeeping on {}.'.format(name))
        data[19:30, 21:37] = 0

    if name == '5BZQJ1619+5256':
        print('Doing a little housekeeping on {}.'.format(name))
        data[42:46, 48:55] = 0
        data[27:32, 35:42] = 0

    if name == '5BZBJ1037+5711':
        print('Doing a little housekeeping on {}.'.format(name))
        data[21:26, 27:32] = 0
        data[55:61, 39:46] = 0

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
    """Fit a two-dimensional Gaussian to the data.
    """
    x, y = xy
    a = ((np.cos(theta) ** 2) / (2 * sigma_x ** 2) +
         (np.sin(theta) ** 2) / (2 * sigma_y ** 2))
    b = (-(np.sin(2 * theta)) / (4 * sigma_x ** 2) +
         (np.sin(2 * theta)) / (4 * sigma_y ** 2))
    c = ((np.sin(theta) ** 2)/(2 * sigma_x ** 2) +
         (np.cos(theta) ** 2) / (2 * sigma_y ** 2))
    g = (offset + amplitude * np.exp(-(a * ((x - x0) ** 2) +
         2 * b * (x - x0) * (y - y0) + c * ((y - y0) ** 2))))
    return g
    # g_ravel = g.ravel()
    # return g_ravel


# def make_model(data, amplitude=1, sigma_x=30, sigma_y=30, theta=0, offset=0):
#     '''Fit a model to the data.'''
#     len_x, len_y = data.shape
#     range_x = range(len_x)
#     range_y = range(len_y)
#     x, y = np.meshgrid(range_x, range_y)
#     x0 = len_x / 2
#     y0 = len_y / 2
#     data_ravel = data.ravel()
#     p0 = (amplitude, x0, y0, sigma_x, sigma_y, theta, offset)
#     popt, pcov = opt.curve_fit(gaussian, (x, y), data_ravel, p0=p0)
#     model = gaussian((x, y), *popt)
#     model_reshape = model.reshape(len_x, len_y)
#     return model_reshape


def match_peaks(data, model, cval=0):
    """Shift the data so that the index of the maximum value in the data
    matches up with that of the model. This is needed to ensure accurate
    subtraction.
    """
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


def jy_per_beam_to_jy(flux, beam=6, sigma=5):
    '''Convert Jansky per beam measurements to Jansky measurements.'''
    bmaj = beam * u.arcsec
    bmin = beam * u.arcsec
    fwhm_to_sigma = 1 / np.sqrt(8 * np.log(2))
    beam_area = 2 * np.pi * bmaj * bmin * (fwhm_to_sigma ** 2)
    beam_area = beam_area.to(u.sr)
    flux_jy_per_beam = flux * u.Jy / u.sr
    flux_jy = flux_jy_per_beam * beam_area
    # BUG this is way off but the code looks right
    # see http://docs.astropy.org/en/stable/api/astropy.units.equivalencies.brightness_temperature.html
    # also see https://astronomy.stackexchange.com/a/20391/20014
    # also see my note from 2019-06-20
    return flux_jy #


def diffuse_fraction(df, name, blazar, diffuse, threshold):
    '''Calculate what the fraction of emission is diffuse.'''
    blazar_catalogue = df.loc[name, 'Total_flux']
    blazar_image = blazar[blazar > threshold].sum()
    # print(f'There is {blazar_image} Jy per beam above the given \N{GREEK SMALL LETTER SIGMA}.')
    diffuse_image = diffuse[diffuse > threshold].sum()
    fraction = diffuse_image / blazar_image
    blazar_diffuse_catalogue = blazar_catalogue * fraction
    blazar_core_catalogue = blazar_catalogue * (1 - fraction)
    print('{0:.2%} of the emission is not from the core.'.format(fraction))
    results = {'total': blazar_catalogue * 1000, 'core': blazar_core_catalogue * 1000, 'diffuse': blazar_diffuse_catalogue * 1000}  # mJy
    print('The total flux density is {total:.2f} mJy, with the core contributing {core:.2f} mJy and the diffuse emission contributing {diffuse:.2f} mJy.'.format(**results))
    return name, blazar_diffuse_catalogue


def make_plot(position, data, title, rows=2, columns=7, origin='lower',
              vmin=0, vmax=1, zmin=0, zmax=1, axis='off', elev=24, azim=29,
              linewidth=0.25, pane_colour='white', line_colour='black',
              cmap='magma_r', projection='3d', stretch='arcsinh', levels='',
              plot='', layer='', contour_colours=['blue', 'magenta'], pad=0.05,
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


def new_plots(pos, data, layer=None, plot_type='blazar', vmax=1, levels=[]):
    ax0 = plt.subplot(1, 2, pos)
    ax0.get_xaxis().set_ticks([])
    ax0.get_yaxis().set_ticks([])
    image_size = 60 * 2  # image is 2 * 2 arcminutes
    len_x, len_y = data.shape
    x0, x1 = len_x * (0.65 + 0.125), len_x * 0.9
    y = len_y * 0.1
    plt.text((x0 + x1) / 2, y * (0.7), str(int(image_size * 0.25 * 0.5)) + '"', horizontalalignment='center', verticalalignment='center', fontsize=15, color='w')
    plt.plot([x0, x1], [y, y], 'w-', lw=1.5)

    image = ax0.imshow(data, origin='lower', cmap='viridis', vmax=vmax, vmin=0, norm=DS9Normalize(stretch='arcsinh'))
    circle1 = plt.Circle([8, 8], (x1 - x0) * (6 / 15) / 2, ls='--', fc='None', ec='w')
    ax0.add_artist(circle1)
    divider = make_axes_locatable(ax0)
    cax0 = divider.append_axes('bottom', size='5%', pad=0.05)
    cbar = plt.colorbar(image, cax=cax0, orientation='horizontal')
    cbar.ax.tick_params(labelsize=15)
    cbar.ax.set_xlabel(r'Jy beam$^{-1}$', fontsize=15)
    if plot_type is 'blazar':
        ax0.contour(data, levels=levels, origin='lower', colors='white', alpha=0.67)
    elif plot_type is 'diffuse':
        ax0.contour(layer, levels=levels, origin='lower', colors='white', alpha=0.67)
        ax0.contour(data, levels=levels, origin='lower', colors='magenta')  # in this case, data is the core subtracted array
        # ax0.contourf(data, levels=[levels[0],levels[0]*100], origin='lower', colors='magenta', alpha=0.2)


def main():
    """Fit a Gaussian point spread function to a point source and subtract it
       from a source with diffuse emission.
    """
    catalogue_dir = '/data5/sean/ldr2/catalogues/'
    # catalogue = f'{catalogue_dir}ldr2-point-sources-near-bllacs.csv'
    csv = f'{catalogue_dir}LDR2 and BZCAT 10_ crossmatch - Copy of BL Lacs.csv'
    output = f'{catalogue_dir}../images/core-subtraction/'

    font = 'STIXGeneral'
    math_font = 'cm'
    font_size = 12
    figsize = (17, 8)
    # testing = False
    new_size = 10
    my_blazars, my_diffuse = [], []

    # df_blazars = get_df(csv, format='csv', index='BZCAT name')
    df_blazars = pd.read_csv(csv)
    df_blazars.set_index('BZCAT name', inplace=True)
    # (blazar_names, blazar_positions,
    #  catalogues, images) = get_position(df_blazars, cat_dir=catalogue_dir)
    for i, (blazar_name, blazar_ra, blazar_dec, bmaj, bmin, angle, peak_flux,
            image, rms) in enumerate(zip(df_blazars.index.tolist(),
                                         df_blazars['BZCAT RA'],
                                         df_blazars['BZCAT Dec'],
                                         df_blazars['Point major'],
                                         df_blazars['Point minor'],
                                         df_blazars['Point angle'],
                                         df_blazars['Peak_flux'],
                                         df_blazars['Mosaic_ID'],
                                         df_blazars['Isl_rms'])):
        blazar_position = [blazar_ra, blazar_dec]
        # if testing:
        #     if i != 0:  # do one at a time
        #         sys.exit()
        # if blazar_name != '5BZQJ1437+3519':
        #     continue
        print(f'Analysing {blazar_name} which is in {image} (BL Lac {i + 1} of'
              f' {len(df_blazars.index.tolist())}).')
        # df_cat = get_df(catalogue, format='csv', index='Source_Name')
        # (point_source_id,
        #  point_source_position = nearest_point_source(df_cat,
        #                                               blazar_position)
        # if point_source_id is False:
        #     print('Going to the next iteration.')
        #     continue
        # hdu, wcs = get_fits(filename=image)
        hdu = fits.open(f'/data5/sean/ldr2/mosaics/{image}-mosaic.fits')[0]
        wcs = WCS(hdu.header, naxis=2)
        sky_position = SkyCoord(blazar_position[0], blazar_position[1],
                                unit='deg')
        if blazar_name == '5BZB J1202+4444':
            size = [3, 3] * u.arcmin
        elif blazar_name == '5BZB J1419+5423':
            size = [4, 4] * u.arcmin
        else:
            size = [2, 2] * u.arcmin
        cutout = Cutout2D(np.squeeze(hdu.data), sky_position, size=size,
                          wcs=wcs)
        blazar_data = cutout.data  # 1 pixel = 1.5"
        blazar_regrid = blazar_data
        # blazar_regrid = regrid(blazar_data, new_size=new_size,
        #                        normalise=False)
        # 1 pixel of blazar_regrid = 1.5" / new_size
        # blazar_regrid = blazar_data
        x0, y0 = np.unravel_index(blazar_regrid.argmax(), blazar_regrid.shape)
        # plt.imshow(blazar_regrid, origin='lower')
        # plt.show()
        # blazar_data = housekeeping(blazar_name, blazar_data)
        # point_source_data = get_data(position=point_source_position, hdu=hdu,
        #                              wcs=wcs)
        # NOTE peak and total values change with regridding
        # point_source_regrid = regrid(point_source_data, new_size=new_size,
        #                              normalise=False)
        # model = make_model(point_source_regrid, sigma_x=4, sigma_y=4)
        # point_source_residual = point_source_regrid - model
        # scaled_model = (model * np.max(blazar_regrid) /
        #                 np.max(point_source_regrid))
        x_, y_ = blazar_regrid.shape
        xy = np.meshgrid(np.linspace(0, x_ - 1, x_),
                         np.linspace(0, y_ - 1, y_))
        scaled_model = gaussian(xy=xy,
                                amplitude=np.max(blazar_regrid),  # jansky
                                x0=x0,
                                y0=y0,
                                sigma_x=bmaj / 2 / 1.5,  # * new_size,
                                sigma_y=bmin / 2 / 1.5,  # * new_size,
                                theta=angle,
                                offset=0)
        # plt.imshow(scaled_model, origin='lower')
        # plt.show()
        plt.figure(figsize=figsize)
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['mathtext.fontset'] = 'dejavuserif'
        matplotlib.rcParams['xtick.major.size'] = 10
        matplotlib.rcParams['xtick.major.width'] = 2
        matplotlib.rcParams['xtick.minor.size'] = 5
        matplotlib.rcParams['xtick.minor.width'] = 2
        matplotlib.rcParams['ytick.major.size'] = 10
        matplotlib.rcParams['ytick.major.width'] = 2
        matplotlib.rcParams['ytick.minor.size'] = 5
        matplotlib.rcParams['ytick.minor.width'] = 2
        matplotlib.rcParams['axes.linewidth'] = 2

        ax0 = plt.subplot(1, 3, 1, projection=wcs)
        ax0.imshow(blazar_regrid, origin='lower', cmap='RdGy',
                   vmax=np.max(blazar_regrid), vmin=-np.max(blazar_regrid))
        # norm=DS9Normalize(stretch='arcsinh'))
        ax0.tick_params(axis='both', which='major', labelsize=20)
        beam = Circle((6, 6), radius=2, linestyle='dashed', lw=2, fc='none',
                      edgecolor='blue')  # radius=2 pixels -> 3" -> diameter=6"
        ax0.add_patch(beam)
        plt.xlabel('Right ascension', fontsize=20, color='black')
        plt.ylabel('Declination', fontsize=20, color='black')

        ax1 = plt.subplot(1, 3, 2, projection=wcs)
        ax1.imshow(scaled_model, origin='lower', cmap='RdGy',
                   vmax=np.max(blazar_regrid), vmin=-np.max(blazar_regrid))
        # norm=DS9Normalize(stretch='arcsinh'))
        ax1.tick_params(axis='both', which='major', labelsize=20)
        beam = Circle((6, 6), radius=2, linestyle='dashed', lw=2, fc='none',
                      edgecolor='blue')  # radius=2 pixels -> 3" -> diameter=6"
        ax1.add_patch(beam)
        plt.xlabel('Right ascension', fontsize=20, color='black')
        plt.ylabel('Declination', fontsize=20, color='black')

        ax2 = plt.subplot(1, 3, 3, projection=wcs)
        ax2.imshow(blazar_regrid - scaled_model, origin='lower',
                   cmap='RdGy', vmin=-np.max(blazar_regrid),
                   vmax=np.max(blazar_regrid))
        # norm=DS9Normalize(stretch='arcsinh'))
        ax2.tick_params(axis='both', which='major', labelsize=20)
        beam = Circle((6, 6), radius=2, linestyle='dashed', lw=2, fc='none',
                      edgecolor='blue')  # radius=2 pixels -> 3" -> diameter=6"
        ax2.add_patch(beam)
        plt.xlabel('Right ascension', fontsize=20, color='black')
        plt.ylabel('Declination', fontsize=20, color='black')

        plt.tight_layout()
        plt.show()
        return
        continue
        if i > 1:
            return
        # blazar_shifted = match_peaks(blazar_regrid, scaled_model)
        # blazar_residual = blazar_shifted - scaled_model
        # blazar_regrid_back = regrid(blazar_shifted, new_size=1 / new_size, normalise=False)  # regrid the blazar and blazar residual data back to the native resolution
        # blazar_residual_regrid_back = regrid(blazar_residual, new_size=1 / new_size, normalise=False)
        # five_sigma = get_noise_catalogue(df_blazars, blazar_name)
        # if blazar_name == '5BZQJ1437+3519':
        #     print('Doing a little more housekeeping on {}.'.format(blazar_name))
        #     blazar_regrid_back[:, 60:] = 0
        #     blazar_residual_regrid_back[:, 60:] = 0
        b, d = diffuse_fraction(df=df_blazars, name=blazar_name, blazar=blazar_regrid_back, diffuse=blazar_residual_regrid_back, threshold=five_sigma)
        my_blazars.append(b)
        my_diffuse.append(d)
        continue  # skip the plotting as I have that already
        savefig = output + '/' + blazar_name + '.png'
        matplotlib.rcParams['font.family'] = font
        matplotlib.rcParams['mathtext.fontset'] = math_font
        matplotlib.rcParams['font.size'] = font_size
        plt.figure(figsize=figsize)
        # make_plot(position=1, data=point_source_regrid, title=point_source_id, vmax=np.max(point_source_regrid))
        # make_plot(position=2, data=model, title=point_source_id + ' model', vmax=np.max(point_source_regrid))
        # make_plot(position=3, data=point_source_residual, title=point_source_id + ' residual', vmax=np.max(point_source_regrid))
        # make_plot(position=4, data=blazar_shifted, title=blazar_name, vmax=np.max(blazar_shifted))
        # make_plot(position=5, data=blazar_residual, title=blazar_name + ' diffuse', vmax=np.max(blazar_shifted))
        # make_plot(position=6, data=blazar_regrid_back, title=blazar_name, levels=five_sigma, plot='blazar', vmax=np.max(blazar_regrid_back))
        # make_plot(position=7, data=blazar_residual_regrid_back, title=blazar_name + ' diffuse', levels=five_sigma, plot='diffuse', layer=blazar_regrid_back, vmax=np.max(blazar_regrid_back))
        # if testing:
        #     plt.show()
        # else:
        #     plt.savefig(savefig, bbox_inches=bbox_inches)
        #     print(f'Done! The plot is saved. View it with this: gpicview {savefig}')
        # TODO make a df from the results and export it as a csv
        # TODO make a new plotter only showing the 2d blazar_regrid_back and blazar_residual_regrid_back with viridis
        #      could I get this contour onto the radio image with proper axes?
        #      my results might not be correct as I am manually blocking out areas for certain sources
        # TODO plot negative flux to see if there are any bowls, using a diverging colour scale
        # TODO does the flux within 5 sigma equal the catalogue flux?
        new_plots(pos=1, data=blazar_regrid_back, vmax=np.max(blazar_regrid_back), levels=five_sigma)
        # old_max = np.max(blazar_regrid_back)
        # blazar_regrid_back[24:43, 36:42] = 0
        # blazar_regrid_back[38:40, 42:43] = 0
        # blazar_regrid_back[35:38, 35:36] = 0
        # blazar_regrid_back[27:30, 42:43] = 0
        # new_plots(pos=2, data=blazar_regrid_back, vmax=old_max, levels=five_sigma, layer=blazar_regrid_back, plot_type='diffuse')
        new_plots(pos=2, data=blazar_residual_regrid_back, vmax=np.max(blazar_regrid_back), levels=five_sigma, layer=blazar_regrid_back, plot_type='diffuse')
        # plt.title(blazar_name, color='white')
        plt.tight_layout()
        plt.show()
        # plt.savefig(f'{output}/2-plots-{blazar_name}.png')

    print()
    for b, d in zip(my_blazars, my_diffuse):
        print(f'{b} {d}')


if __name__ == '__main__':
    main()
