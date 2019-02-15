#!/usr/bin/env python3.6

'''Plot one FITS file (e.g. LOFAR radio data) over another (e.g. SDSS optical
data).'''

import warnings
warnings.filterwarnings('ignore')  # supress warnings

import aplpy
import argparse
import matplotlib.pyplot as plt
import montage_wrapper as montage
import numpy as np
import os
from PIL import Image, ImageDraw

__author__ = 'Sean Mooney'
__date__ = '13 February 2019'

def get_image(fits, centre, cmap, vmin, vmax, radius=450, dpi=256, max_dpi=0,
              format='png'):
    '''Get a png image from a FITS file. The size is the width of the square in
    arcseconds.'''

    plt.rcParams['axes.edgecolor'] = 'white'
    image = aplpy.FITSFigure(fits)  # north=True
    image.recenter(centre[0], centre[1], radius=radius / 60 / 60)
    image.show_colorscale(cmap=cmap, vmin=vmin, vmax=vmax)
    image.frame.set_linewidth(0)
    image.hide_axis_labels()
    image.hide_tick_labels()

    save = os.path.splitext(fits)[0] + '.' + format
    image.save(save, dpi=dpi, max_dpi=max_dpi)

    return save


def white_to_transparent(image):
    '''Makes white pixels of an image transparent.'''

    x = np.asarray(image.convert('RGBA')).copy()

    x[:, :, 3] = (255 * (x[:, :, :3] != 255).any(axis=2)).astype(np.uint8)

    return Image.fromarray(x)


def overlay_image(front, back, size=1024, alpha=0.5):
    '''Takes two images and plots the front one with moderate transparency over
    the back image.'''

    image1 = Image.open(front)
    image2 = Image.open(back)
    transparent1 = image1 # white_to_transparent(image1)  # (0, 0, 0, 0) is white and transparent
    output = Image.blend(transparent1, image2, alpha=alpha)  # alpha=0 shows transparent1, alpha=1 shows image2

    save = '/'.join(front.split('/')[:-1]) + '/output.png'
    output.save(save)
    print('The output image is at {}.'.format(save))

    os.remove(front)  # delete intermediate files
    os.remove(back)


def main():
    '''Plot one FITS file (e.g. LOFAR radio data) over another (e.g. SDSS
    optical data).'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-f', '--radio', required=False, type=str,
                        default='/home/sean/Downloads/fits/m51-r.fits',
                        help='First FITS file to be plotted.')

    parser.add_argument('-F', '--optical', required=False, type=str,
                        default='/home/sean/Downloads/fits/m51-o.fits',
                        help='Second FITS file to be plotted.')

    parser.add_argument('-r', '--ra', required=False, type=float,
                        default=202.484167,
                        help='RA of the source.')

    parser.add_argument('-d', '--dec', required=False, type=float,
                        default=47.230556,
                        help='Declination of the source.')

    args = parser.parse_args()
    radio = args.radio
    optical = args.optical
    ra = args.ra
    dec = args.dec

    radio_image = get_image(fits=radio, centre=[ra, dec], cmap='hot', vmin=0.001, vmax=0.005)
    optical_image = get_image(fits=optical, centre=[ra, dec], cmap='gray', vmin=3000, vmax=10000)
    overlay_image(front=radio_image, back=optical_image)


if __name__ == '__main__':
    main()
