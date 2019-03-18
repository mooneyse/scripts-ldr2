#!/usr/bin/env python3

'''Fit a model of non-thermal radiation from a population of relativistic
particles to an observed spectrum.'''

import argparse
import astropy.units as u
from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
from naima.models import PowerLaw, Synchrotron

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '06 March 2019'


def model(csv, B=100 * u.uG):
    '''Build the synchrotron model.'''

    # the particle distribution function takes electron energies as an array or
    # float and returns the particle energy density in units of number of
    # electrons per unit energy as an array or float
    powerlaw = PowerLaw(amplitude=1e36 / u.eV, e_0=1 * u.TeV, alpha=1)
    synchrotron = Synchrotron(particle_distribution=powerlaw, B=B)

    # data = ascii.read(csv)
    # figure = naima.plot_data(data, e_unit=u.eV)

    # plot the computed model emission
    fig, ax = plt.subplots()
    energy = np.logspace(-7, 15, 100) * u.eV
    ax.loglog(energy, synchrotron.sed(energy, 1 * u.kpc))
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel(r'$E^2 dN/dE$')
    ax.set_xlim(1e-7, 1e9)
    ax.set_ylim(1e-16, 1e0)
    plt.show()


def main():
    '''Fit a model of non-thermal radiation from a population of relativistic
    particles to an observed spectrum.'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-c', '--csv', type=str, help='CSV containing data',
                        default='/home/sean/Downloads/workbooks/oj287.csv')

    args = parser.parse_args()
    csv = args.csv

    model(csv)


if __name__ == '__main__':
    main()
