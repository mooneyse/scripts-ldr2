#!/usr/bin/env python3

'''Fit a Gaussian point spread function to a point source and subtract it from
a source with diffuse emission.'''

from numpy import *
from scipy import optimize
from pylab import *
import random

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: ravel(gaussian(*p)(*indices(data.shape)) - data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

# Create the gaussian data
Xin, Yin = mgrid[0:201, 0:201]
data = gaussian(3, 100, 100, 20, 40)(Xin, Yin) + random.random()

matshow(data, cmap=cm.gist_earth_r)

params = fitgaussian(data)
print(type(data))
fit = gaussian(*params)

contour(fit(*indices(data.shape)), cmap=cm.copper)
ax = gca()
(height, x, y, width_x, width_y) = params

show()

#--------------------------------------


import argparse
import numpy as np
from matplotlib import rcParams
from astropy.table import Table
from photutils.datasets import make_random_gaussians_table, make_noise_image, make_gaussian_sources_image
from photutils.detection import IRAFStarFinder
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
from photutils.psf import IterativelySubtractedPSFPhotometry
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
import matplotlib.pyplot as plt

sigma_psf = 2.0
sources = Table()
sources['flux'] = [700, 800, 700, 800]
sources['x_mean'] = [12, 17, 12, 17]
sources['y_mean'] = [15, 15, 20, 20]
sources['x_stddev'] = sigma_psf*np.ones(4)
sources['y_stddev'] = sources['x_stddev']
sources['theta'] = [0, 0, 0, 0]
sources['id'] = [1, 2, 3, 4]
tshape = (32, 32)
image = (make_gaussian_sources_image(tshape, sources) +
         make_noise_image(tshape, type='poisson', mean=6.,
                          random_state=1) +
         make_noise_image(tshape, type='gaussian', mean=0.,
                          stddev=2., random_state=1))

bkgrms = MADStdBackgroundRMS()
std = bkgrms(image)
iraffind = IRAFStarFinder(threshold=3.5*std,
                          fwhm=sigma_psf*gaussian_sigma_to_fwhm,
                          minsep_fwhm=0.01, roundhi=5.0, roundlo=-5.0,
                          sharplo=0.0, sharphi=2.0)
daogroup = DAOGroup(2.0*sigma_psf*gaussian_sigma_to_fwhm)
mmm_bkg = MMMBackground()
fitter = LevMarLSQFitter()
psf_model = IntegratedGaussianPRF(sigma=sigma_psf)
photometry = IterativelySubtractedPSFPhotometry(finder=iraffind,
                                                group_maker=daogroup,
                                                bkg_estimator=mmm_bkg,
                                                psf_model=psf_model,
                                                fitter=LevMarLSQFitter(),
                                                niters=1, fitshape=(11,11))
result_tab = photometry(image=image)
residual_image = photometry.get_residual_image()

plt.subplot(1, 2, 1)
plt.imshow(image, cmap='viridis', aspect=1, interpolation='nearest',
           origin='lower', vmin=0, vmax=40)
plt.title('Simulated data')
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
plt.subplot(1 ,2, 2)
plt.imshow(residual_image, cmap='viridis', aspect=1,
           interpolation='nearest', origin='lower', vmin=0, vmax=40)
plt.title('Residual Image')
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
plt.show()




__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '20 March 2019'

def subtract_gaussian():
    '''Subtract Gaussian.'''

    pass


def main():
    '''Fit a Gaussian point spread function to a point source and subtract it
    from a source with diffuse emission.'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-p',
                        '--point',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/deep-fields/point.fits',
                        help='FITS file containing a point source')

    parser.add_argument('-d',
                        '--diffuse',
                        required=False,
                        type=str,
                        default='/mnt/closet/ldr2-blazars/deep-fields/diffuse.fits',
                        help='FITS file containing a diffuse source')

    args = parser.parse_args()
    point = args.point
    diffuse = args.diffuse

    subtract_gaussian()


if __name__ == '__main__':
    main()
