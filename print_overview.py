#!/usr/bin/env python3

import numpy as np
import pandas as pd
from uncertainties import ufloat


def get_mean(dataset, column, errors):
    quantity = []
    for value, error in zip(dataset[column], dataset[errors]):
        if value == '-':  # usually '-'
            quantity.append(ufloat(np.nan, np.nan))
        else:
            if type(value) is not float:  # core dominance read in a as string
                value = float(value)
            if type(error) is not float:
                error = float(error)
            quantity.append(ufloat(value, error))
    average_quantity = np.sum(quantity) / len(quantity)
    return f'{average_quantity.n} ± {average_quantity.s}'


df = pd.read_csv('LDR2 and BZCAT 10_ crossmatch - Basic.csv')

extended = df[df['Compact'] == False]
unresolved = df[df['Compact'] == True]
gamma = df[df['4LAC name'] != '-']
non_gamma = df[df['4LAC name'] == '-']

extended_gamma = df[(df['Compact'] == False) & (df['4LAC name'] != '-')]
extended_non_gamma = df[(df['Compact'] == False) & (df['4LAC name'] == '-')]

unresolved_gamma =df[(df['Compact'] == True) & (df['4LAC name'] != '-')]
unresolved_non_gamma = df[(df['Compact'] == True) & (df['4LAC name'] == '-')]

my_dict = {'total': df,
           'extended': extended,
           'unresolved': unresolved,
           'γ': gamma,
           'non-γ': non_gamma,
           'extended γ': extended_gamma,
           'extended non-γ': extended_non_gamma,
           'unresolved γ': unresolved_gamma,
           'unresolved non-γ': unresolved_non_gamma}

for label, dataset in my_dict.items():
    results = {}
    results['N'] = len(dataset.index)
    results['z'] = np.average(dataset['redshift'])
    results['L'] = get_mean(dataset,
                            column='Luminosity with FIRST index (W/Hz)',
                            errors='Luminosity error with FIRST index (W/Hz)')
    results['α'] = get_mean(dataset,
                            column='LDR2-to-FIRST index',
                            errors='LDR2-to-FIRST index error')
    results['log(R)'] = get_mean(dataset,
                                 column='Log10(core dominance)',
                                 errors='Log10(core dominance) error')
    results['E'] = get_mean(dataset,
                            column='Extent (kpc)',
                            errors='Extent error (kpc)')
    results['S_γ'] = get_mean(dataset,
                              column='Flux1000',
                              errors='Unc_Flux1000')

    print(f'\n{label}')
    for quantity, value in results.items():
        print(f'{quantity}\t=\t{value}')
