#!/usr/bin/env python3

import os
import numpy as np
from astropy.io import fits
from astropy.stats import mad_std
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse
from scipy.optimize import curve_fit

fits_file = '/mnt/closet/ldr2-blazars/deep-fields/bootes-image.fits'

data = np.squeeze(fits.open(fits_file)[0].data)

field = data[13543:15543, 11523:13523]
plt.subplot(1, 3, 1)
plt.imshow(field, cmap='viridis', vmin=0, vmax=0.01, origin='lower')
plt.title('Blazar field')
plt.axis('off')

diffuse = field[986:1018, 983:1015]
diffuse = diffuse / np.max(diffuse)  # normalise
plt.subplot(1, 3, 2)
plt.imshow(diffuse, cmap='viridis', vmin=0, vmax=1, origin='lower')
plt.title('Blazar source')
plt.axis('off')

point = field[240:272, 336:368]
point = point / np.max(point)  # normalise
plt.subplot(1, 3, 3)
plt.imshow(point, cmap='viridis', vmin=0, vmax=1, origin='lower')
plt.title('Point source')
plt.axis('off')

plt.show()

def gaussian(xy, x0, y0, sigma, amplitude):
    x, y = xy
    A = 1 / (2 * sigma ** 2)
    return  amplitude * np.exp(-A * ((x - x0) ** 2 + (y - y0) ** 2))


def fit(image):
    med = np.median(image)
    image = image - med
    max_index = np.where(image >= np.max(image))
    x0 = max_index[1] # middle of x axis
    y0 = max_index[0] # middle of y axis
    x = np.arange(0, image.shape[1], 1)
    y = np.arange(0, image.shape[0], 1)
    xx, yy = np.meshgrid(x, y) # creates a grid to plot the function over
    sigma = np.std(image)
    amp = np.max(image)
    guess = [x0, y0, sigma, amp]
    lower = [0, 0, 0, 0] #start of data array
    upper = [image.shape[0], image.shape[1], np.max(image), np.max(image) * 2]
    bounds = [lower, upper]
    params, pcov = curve_fit(gaussian, (xx.ravel(), yy.ravel()), image.ravel(), p0=guess, bounds=bounds)
    print(pcov)
    return params


def plotting(image, params):
    fig, ax = plt.subplots()
    ax.imshow(image)
    ax.scatter(params[0], params[1],s = 10, c = 'red', marker = 'x')
    # circle = Circle((params[0], params[1]), 10 * params[2], facecolor = 'none', edgecolor = 'white', linewidth = 1)
    circle = Ellipse((params[0], params[1]), 10 * params[2], 10 * params[3], facecolor = 'none', edgecolor = 'white', linewidth = 1)

    ax.add_patch(circle)
    plt.show()

parameters = fit(point)
plotting(point, parameters)






# from matplotlib import cm
# import numpy
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# # Set up grid and test data
# nx, ny = diffuse.shape
# x = range(nx)
# y = range(ny)
# data = diffuse  # numpy.random.random((nx, ny))
# hf = plt.figure()
# ha = hf.add_subplot(111, projection='3d')
# X, Y = numpy.meshgrid(x, y)
# ha.plot_surface(X, Y, data, rstride=3, cstride=3, linewidth=1, antialiased=True,
#                 cmap=cm.gray)
# cset = ha.contourf(X, Y, data, zdir='z', offset=-0.5, cmap=cm.viridis)
# ha.set_zlim(-0.5,1)
# ha.set_zticks(np.linspace(0,1,5))
# ha.view_init(27, -21)
# plt.show()
