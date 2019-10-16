#!/usr/bin/env python3

"""Plot PyBDSF componets for each source."""

import matplotlib
matplotlib.use('Agg')
# import os
# import sys
import matplotlib.pylab as pl
import pandas as pd
import aplpy

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '18 June 2019'


def get_ellipses(my_dir):
    """Get the ellipse parameters from the CSV.

    Parameters
    ----------
    csv : string
        Filepath of the CSV containing the Gaussians fitted by PyBDSF.

    Returns
    -------
    tuple
        Associated source name, right ascension, declination, semi-major axis,
        semi-minor axis, position angle, total flux density (Jy) and fractional
        flux density (Jy) for each row in the CSV.
    """
    csv = my_dir + 'catalogues/deep.fields.29.07.2019.gaul.csv'
    df = pd.read_csv(csv)
    df = df.sort_values('Total_flux_2', ascending=False)

    source_name = df['Source name']
    ra = df['RA_2']
    dec = df['DEC_2']
    major = df['Maj_2']  # Maj_2, Maj_img_plane_2, DC_Maj_2, DC_Maj_img_plane_2
    minor = df['Min_2']  # Min_2, Min_img_plane_2, DC_Min_2, DC_Min_img_plane_2
    pa = df['PA_2']  # PA_2, PA_img_plane_2, DC_PA_2, DC_PA_img_plane_2
    flux = df['Total_flux_2']
    fraction = flux / df['Total_flux_1']

    ellipses = []
    for i, _ in enumerate(source_name):
        ellipse = (source_name[i], flux[i], fraction[i],
                   ra[i], dec[i], major[i], minor[i], pa[i])
        ellipses.append(ellipse)

    return ellipses


def get_info(my_dir):
    """Get the unique list of blazar names.

    Parameters
    ----------
    ellipses : list
        Tuples with the ellipse parameters.

    Returns
    -------
    list
        Names of the blazars.
    """
    csv = my_dir + 'catalogues/deep.fields.29.07.2019.gaul.csv'
    df = pd.read_csv(csv)

    names = list(dict.fromkeys(df['Source name']).keys())
    ras = list(dict.fromkeys(df['RA_1']).keys())
    decs = list(dict.fromkeys(df['DEC_1']).keys())
    fluxes = list(dict.fromkeys(df['Peak_flux_1']).keys())

    return names, ras, decs, fluxes


def unpack(s):
    """Convenience function to get a list as a string without the braces.

    Parameters
    ----------
    s : list
        The list to be turned into a string.

    Returns
    -------
    string
        The list as a string.
    """
    x = " ".join(map(str, s))
    return x


def make_region_files(names, ellipses, my_dir):
    """Create region files for each source.

    Parameters
    ----------
    names : list
        Blazar names.
    ellipses : list
        Tuples of the ellipse parameters.

    Returns
    -------
    list
        Names of the created region files.
    """
    header = ('# Region file format: DS9 version 4.1\nglobal color=green das' +
              'hlist=8 3 width=1 font="helvetica 10 normal roman" select=1 h' +
              'ighlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 sou' +
              'rce=1\nicrs\n')
    my_dir = my_dir + 'images/regions/'

    region_files = []
    for name in names:
        region = f'{my_dir}{name}.reg'
        region_files.append(region)
        with open(region, 'w+') as the_file:
            the_file.write(header)

        for name_, flux, fraction, *ellipse in ellipses:
            if name == name_:
                ellipse = f'ellipse {unpack(ellipse)} # color=white width=2\n'
                with open(region, 'a') as the_file:
                    the_file.write(ellipse)

    return region_files


def fits_file(name, ra, my_dir):
    """Get the FITS file for a given source.

    Parameters
    ----------
    name : string
        Source name.
    ra : float
        Right ascension.
    my_dir : string
        Working directory.

    Returns
    -------
    string
        Name of the FITS file of the field the source is in.
    """

    if ra > 156 and ra < 168:
        field_file = f'{my_dir}catalogues/lockman.hole.11.06.2019.img.fits'
    elif ra > 214 and ra < 222:
        field_file = f'{my_dir}catalogues/bootes.11.06.2019.img.fits'
    elif ra > 237 and ra < 248:
        field_file = f'{my_dir}catalogues/elias.n1.11.06.2019.img.fits'
    else:
        raise ValueError(f'{name} not in any field.')
    return field_file


def plot_ellipses(names, ras, decs, peaks, my_dir, ellipses):
    """Use AplPy's show_regions to plot the region files.

    Parameters
    ----------
    names : list
        Source names.
    ras : list
        Right ascensions.
    decs : list
        Declinations.
    peaks : list
        Peak fluxes.
    my_dir : string
        Working directory.
    ellipses : list
        Tuple of ellipse parameters.

    Returns
    -------
    None.
    """
    for name, ra, dec, peak in zip(names, ras, decs, peaks):
        print(f'Imaging {name}.')
        field_file = fits_file(name=name, ra=ra, my_dir=my_dir)
        image = aplpy.FITSFigure(field_file)
        image.recenter(ra, dec, radius=1 / 60)
        for n, _, flux, ra, dec, maj, min, pa in ellipses:
            if name == n:
                col = pl.cm.autumn(int(255 * flux))
                facecolor = tuple((col[0], col[1], col[2], 0.25))
                edgecolor = tuple((col[0], col[1], col[2], 0.5))
                image.show_ellipses(ra, dec, maj * 2, min * 2, pa + 90,
                                    facecolor=facecolor,
                                    edgecolor=edgecolor)

        image.show_colorscale(cmap='viridis', vmin=0, vmax=peak,
                              stretch='arcsinh')
        save = f'{my_dir}images/component-summation/{name}.png'
        image.save(save)
        print(f'{save} done.')
        # command = f'gpicview {my_dir}images/component-summation/{name}.png'
        # os.system(command)
        # sys.exit()


def do_the_math():
    """Find the core and diffuse flux density.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    pass


def main():
    """Plot PyBDSF componets for each source.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    server = True
    my_dir = 'data5/sean' if server else 'mnt/closet'
    my_dir = f'/{my_dir}/deep-fields/'

    ellipses = get_ellipses(my_dir=my_dir)
    names, ras, decs, peaks = get_info(my_dir=my_dir)

    # region_files = make_region_files(namesr=names, ellipses=ellipses,
    #                                  my_dir=my_dir)

    plot_ellipses(names=names, ras=ras, decs=decs, peaks=peaks, my_dir=my_dir,
                  ellipses=ellipses)

    # NOTE could make it so it plots the smallest flux first, by sorting the
    #      ellipses by the fraction flux

    do_the_math()


if __name__ == '__main__':
    main()
