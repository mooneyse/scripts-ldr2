#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.modeling.models import Gaussian2D
from photutils.datasets import make_noise_image
from photutils.isophote import build_ellipse_model, Ellipse, EllipseGeometry
from photutils import EllipticalAperture

fits_file = '/mnt/closet/ldr2-blazars/deep-fields/bootes-image.fits'
my_array = np.squeeze(fits.open(fits_file)[0].data)
field = my_array[13543:15543, 11523:13523]
point = field[240:272, 336:368]
point = my_array[13783:13815, 11859:11891]
data = point / np.max(point)

geometry = EllipseGeometry(x0=16, y0=16, sma=5, eps=0, pa=0)
aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,
                           geometry.sma * (1 - geometry.eps), geometry.pa)

ellipse = Ellipse(data, geometry)
isolist = ellipse.fit_image()

model_image = build_ellipse_model(data.shape, isolist)
residual = data - model_image

fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3)
ax1.imshow(data, origin='lower', vmin=0, vmax=1)
ax1.set_title('Data')
ax1.axis('off')

ax2.imshow(model_image, origin='lower', vmin=0, vmax=1)
ax2.set_title('Model')
ax2.axis('off')

ax3.imshow(residual, origin='lower', vmin=0, vmax=1)
ax3.set_title('Residual')
ax3.axis('off')

plt.show()
