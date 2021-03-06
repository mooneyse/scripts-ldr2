#!/usr/bin/env python3

'''Fit a Gaussian point spread function to a point source and subtract it from
a source with diffuse emission.'''

import warnings
warnings.filterwarnings('ignore')  # supress warnings
from uncertainties import ufloat
from uncertainties.umath import log, sqrt, exp
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
    df_blazars = pd.read_csv(f'{catalogue_dir}LDR2 and BZCAT 10_ crossmatch - '
                             'Copy of BL Lacs.csv')
    df_blazars.set_index('BZCAT name', inplace=True)
    print(f'Index (of {len(df_blazars.index.tolist())}), Source, Total (mJy), '
          'Core (mJy), Diffuse (mJy), Unresolved?')

    for i, (blazar_name, blazar_ra, blazar_dec, bmaj, bmin, angle, peak_flux,
            image, rms, total_flux,
            unresolved) in enumerate(zip(df_blazars['Name'],
                                         df_blazars['BZCAT RA'],
                                         df_blazars['BZCAT Dec'],
                                         df_blazars['Point major'],
                                         df_blazars['Point minor'],
                                         df_blazars['Point angle'],
                                         df_blazars['Peak_flux'],
                                         df_blazars['Mosaic_ID'],
                                         df_blazars['Isl_rms'],
                                         df_blazars['Total_flux'],
                                         df_blazars['Compact'])):
        blazar_position = [blazar_ra, blazar_dec]
        hdu = fits.open(f'/data5/sean/ldr2/mosaics/{image}-mosaic.fits')[0]
        wcs = WCS(hdu.header, naxis=2)
        sky_position = SkyCoord(blazar_position[0], blazar_position[1],
                                unit='deg')
        if blazar_name == 'J1202+4444' or blazar_name == '5BZBJ1325+4115':
            size = [3, 3] * u.arcmin
        elif (blazar_name == 'J1419+5423' or
              blazar_name == '5BZBJ0945+5757'):
            size = [4, 4] * u.arcmin
        else:
            size = [2, 2] * u.arcmin

        cutout = Cutout2D(np.squeeze(hdu.data), sky_position, size=size,
                          wcs=wcs)
        blazar_data = cutout.data  # 1 pixel = 1.5"
        blazar_regrid = blazar_data
        x0, y0 = np.unravel_index(blazar_regrid.argmax(), blazar_regrid.shape)
        x_, y_ = blazar_regrid.shape
        xy = np.meshgrid(np.linspace(0, x_ - 1, x_),
                         np.linspace(0, y_ - 1, y_))
        scaled_model = gaussian(xy=xy,
                                amplitude=peak_flux / 1000,
                                # amplitude=np.max(blazar_regrid),
                                x0=y0,
                                y0=x0,
                                # fwhm from asec to pixels to sigma
                                sigma_x=bmaj / 1.5 / 2.355,
                                sigma_y=bmin / 1.5 / 2.355,
                                theta=angle + 90,
                                offset=0)

        core_flux = np.sum(scaled_model) * 1000 / 18.1294  # divide beam area
        # log_core_dominance = np.log10()
        print(f'{i + 1}, {blazar_name}, {total_flux}, {core_flux}, '
              f'{total_flux - core_flux}, {unresolved}')  #, {log_core_dominance}')

        plt.figure(figsize=(32, 8))
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

        p = 'unresolved' if unresolved else 'resolved'
        plt.suptitle(f'{blazar_name} ({p})', fontsize=30)

        ax0 = plt.subplot(1, 3, 1, projection=wcs)
        ax0.imshow(blazar_regrid, origin='lower', cmap='RdGy',
                   vmax=np.max(blazar_regrid), vmin=-np.max(blazar_regrid))
        # norm=DS9Normalize(stretch='arcsinh'))
        ax0.tick_params(axis='both', which='major', labelsize=20,
                        direction='in')
        beam = Circle((6, 6), radius=2, linestyle='dashed', lw=2, fc='none',
                      edgecolor='black')
        ax0.add_patch(beam)
        plt.xlabel('Right ascension', fontsize=20, color='black')
        plt.ylabel('Declination', fontsize=20, color='black')
        plt.title(r'$S_\mathrm{int}$' + f' = {total_flux:.0f} mJy',
                  fontsize=20)
        ax0.contour(blazar_regrid, levels=[rms * 4 / 1000, rms * 8 / 1000,
                                           rms * 16 / 1000, rms * 32 / 1000],
                    origin='lower', colors=['r', 'magenta', 'yellow', 'cyan'])

        ax1 = plt.subplot(1, 3, 2, projection=wcs)
        ax1.imshow(scaled_model, origin='lower', cmap='RdGy',
                   vmax=np.max(blazar_regrid), vmin=-np.max(blazar_regrid))
        ax1.tick_params(axis='both', which='major', labelsize=20,
                        direction='in')
        beam = Circle((6, 6), radius=2, linestyle='dashed', lw=2, fc='none',
                      edgecolor='black')
        ax1.add_patch(beam)
        plt.xlabel('Right ascension', fontsize=20, color='black')
        plt.ylabel('Declination', fontsize=20, color='black')
        plt.title(r'$S_\mathrm{core}$' + f' = {core_flux:.0f} mJy',
                  fontsize=20)
        ax1.contour(blazar_regrid, levels=[rms * 4 / 1000], origin='lower',
                    colors='r')
        ax1.contour(scaled_model, levels=[rms * 4 / 1000], origin='lower',
                    colors='g')

        ax2 = plt.subplot(1, 3, 3, projection=wcs)
        ax2.imshow(blazar_regrid - scaled_model, origin='lower',
                   cmap='RdGy', vmin=-np.max(blazar_regrid),
                   vmax=np.max(blazar_regrid))
        ax2.tick_params(axis='both', which='major', labelsize=20,
                        direction='in')
        beam = Circle((6, 6), radius=2, linestyle='dashed', lw=2, fc='none',
                      edgecolor='black')
        ax2.add_patch(beam)
        plt.xlabel('Right ascension', fontsize=20, color='black')
        plt.ylabel('Declination', fontsize=20, color='black')
        plt.title(r'$S_\mathrm{diffuse}$' + f' = {total_flux - core_flux:.0f}'
                  + ' mJy', fontsize=20)
        ax2.contour(blazar_regrid, levels=[rms * 4 / 1000], origin='lower',
                    colors='r')
        ax2.contour(blazar_regrid - scaled_model, levels=[rms * 4 / 1000],
                    origin='lower', colors='blue', ls='dashed')

        # if blazar_name == '5BZB J1340+4410':
        #     plt.show()
        #     print(blazar_name, 'blazar peak:', x0, y0)
        #     x0, y0 = np.unravel_index(scaled_model.argmax(),
        #                               scaled_model.shape)
        #     print(blazar_name, 'Gaussian peak:', x0, y0)
        # else:
        plt.savefig(f'{catalogue_dir}../images/core-subtraction'
                    f'/core-sub-{blazar_name.replace(" ", "")}.png')
        plt.close()


if __name__ == '__main__':
    main()
