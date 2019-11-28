#!/usr/bin/env python3

from uncertainties import ufloat
from uncertainties.umath import log, sqrt, exp, log10
import pandas as pd
import numpy as np

print('Name, Core dominance, Core dominance error, log10(core dom), log10(core'
      ' dom error)')
df = pd.read_csv('LDR2 and BZCAT 10_ crossmatch - Basic.csv')

for name, core, total, tot_err, unresolved in zip(df['Name'],
                                                  df['Raw core calculation'],
                                                  df['Total_flux'],
                                                  df['E_Total_flux'],
                                                  df['Compact']):
    if unresolved:
        print(f'{name}, -, -, -, -')
        continue
    S_core = ufloat(core, core * (tot_err / total))
    S_int = ufloat(total, tot_err)
    core_dom = S_core / (S_int - S_core)
    log_core_dom = log10(core_dom)
    # print(S_core, S_int, core_dom)
    print(f'{name}, {core_dom.n}, {core_dom.s}, {log_core_dom.n},'
          f'{log_core_dom.s}')


def alpha(flux1, flux2, freq1=ufloat(144e6, 24e6),
          freq2=ufloat(1400e6, 0)):
    """Get spectral index with error. Default is LDR2 and FIRST.
    """
    return log(flux1 / flux2) / log(freq1 / freq2)


def get_dl_and_kpc_per_asec(z, H0=70, WM=0.26, WV=0.74):
    """Get luminosity distance. See
    """
    WR = 0.  # Omega(radiation)
    WK = 0.  # Omega curvaturve = 1-Omega(total)
    c = 299792.458  # velocity of light in km/sec
    DTT = 0.5  # time from z to now in units of 1/H0
    age = 0.5  # age of Universe in units of 1/H0
    zage = 0.1  # age of Universe at redshift z in units of 1/H0
    DCMR = 0.0  # comoving radial distance in units of c/H0
    DA = 0.0  # angular size distance
    DA_Mpc = 0.0
    kpc_DA = 0.0
    DL = 0.0  # luminosity distance
    DL_Mpc = 0.0
    a = 1.0  # 1/(1+z), the scale factor of the Universe
    az = 0.5  # 1/(1+z(object))
    h = H0 / 100.
    WR = 4.165E-5 / (h * h)  # includes 3 massless neutrino species,
    WK = 1 - WM - WR - WV
    az = 1.0 / (1 + 1.0 * z)
    age = 0.
    n = 1000  # number of points in integrals
    for i in range(n):
        a = az * (i + 0.5) / n
        adot = sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
        age = age + 1. / adot
    zage = az * age / n
    DTT = 0.0
    DCMR = 0.0
    for i in range(n):
        a = az + (1 - az) * (i + 0.5) / n
        adot = sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
        DTT = DTT + 1. / adot
        DCMR = DCMR + 1. / (a * adot)
    DTT = (1. - az) * DTT / n
    DCMR = (1. - az) * DCMR / n
    age = DTT + zage
    ratio = 1.00
    x = sqrt(abs(WK)) * DCMR
    if x > 0.1:
        if WK > 0:
            ratio = 0.5 * (exp(x) - exp(-x)) / x
        else:
            ratio = sin(x) / x
    else:
        y = x * x
        if WK < 0:
            y = -y
        ratio = 1. + y / 6. + y * y / 120.
    return (c / H0) * ((az * ratio * DCMR) / (az * az)) * 3.086e22, ((c / H0) * az * ratio * DCMR) / 206.264806


def luminosity(S_obs, z, D_L=0, alpha=0):
    """Get radio luminosity with error. Default is LDR2. See
    https://www.fxsolver.com/browse/formulas/Radio+luminosity.
    """
    if D_L == 0:
        D_L, _ = get_dl_and_kpc_per_asec(z=z)
    return (S_obs * 4 * np.pi * (D_L ** 2)) / (1 + z) ** ( 1 + alpha)


df = pd.read_csv('/home/sean/Downloads/LDR2 and BZCAT 10_ crossmatch - Sheet7.csv')
spindices = {}
i1, i2 = [], []
print('\nSpectral index:')

for source, lotss, lotss_err, first, first_err, nvss in zip(df['Name'],
                                                            df['Total_flux'],
                                                            df['E_Total_flux'],
                                                            df['FINT'],
                                                            df['RMS'],
                                                            df['S1.4.1']):
    spindex1 = alpha(flux1=ufloat(lotss, lotss * 0.2),  # lotss_err
                     flux2=ufloat(first, first * 0.1))
    spindex2 = alpha(flux1=ufloat(lotss, lotss * 0.2),  # lotss_err
                     flux2=ufloat(nvss, nvss * 0.1))
#    spindices[source] = spindex
    print(f'{source},{spindex1.n},{spindex1.s},{spindex2.n},{spindex2.s}')
#    print(f'{source},{lotss},{first},{nvss}')
#    print(f'{source}: {spindex}')
    i1.append(spindex1)
    i2.append(spindex2)

#print('\nLuminosity:')
#print(df['S1.4.1'])
for source, flux, flux_err, z, alpha in zip(df['Name'], df['Total_flux'],
                                            df['E_Total_flux'], df['z'],
                                            spindices):
    power = luminosity(S_obs=ufloat(flux * 1e-29, flux * 0.2 * 1e-29), z=z,
            alpha=spindices[source])

#    print(f'{source}: {power}')

print('Averages:')
print('FIRST:', sum(i1)/len(i1))
print('NVSS:', sum(i2)/len(i2))

print('\nAverage flux density:')
lotss, first, nvss = [], [], []
for l, le, f, fe, n, ne in zip(df['Total_flux'], df['E_Total_flux'],
                               df['FINT'], df['RMS'], df['S1.4.1'],
                               df['e_S1.4']):
    lotss.append(ufloat(l, l*0.2))
    first.append(ufloat(f, f*0.1))
    nvss.append(ufloat(n, n*0.1))

print('LDR2:', sum(lotss)/len(lotss))
print('FIRST:', sum(first)/len(first))
print('NVSS:', sum(nvss)/len(nvss))
