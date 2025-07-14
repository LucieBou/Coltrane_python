# -*- coding: utf-8 -*-

'''
Coltrane - timing_combinations

@author Lucie Bourreau & Neil Banas
@date 2023/08/21
'''

import numpy as np

def timing_combinations(forcing, p):
    '''
    Figures out the full set of values of t0 (spawning dates) and the fields of s (strategy) to consider.
    
    Parameters
    ----------
    forcing: dict
        Set of forcing, composed of prey and temperature (surface & deep) cycle over several years.
    p: dict
        Set of parameters.
        
    Returns
    -------
    t0: array
        Spawning dates.
    s: dict
        Timing strategy vector: diapause exit date, diapause entry date and date of egg production.
    '''
    
    ## SPAWNING DATE
    t0 = np.arange(forcing['t'][0], forcing['t'][-1]-365 + .1 + p['dt_spawn'], p['dt_spawn'])
    
    ## YEARDAY OF DIAPAUSE EXIT
    tdia_exit = p['tdia_exit']
    if isinstance(tdia_exit, list):
        if len(tdia_exit) == 0 or (tdia_exit == 'default'):
            tdia_exit = np.arange(0, 365 / 2, p['dt_dia'])

    ## YEARDAY OF DIAPAUSE ENTRY
    tdia_enter = p['tdia_enter']
    if isinstance(tdia_enter, list):
        if len(tdia_enter) == 0 or tdia_enter == 'default':
            tdia_enter = np.arange(max(tdia_exit) + p['dt_dia'], 365, p['dt_dia'])
        
    ## THE DATE THAT EGG PRODUCTION BEGINS RELATIVE TO t0
    dtegg = p['dtegg']
    if len(dtegg) == 0:
        dteggmin = (p['min_genlength_years'] - 0.5) * 365
        dteggmin = max(dteggmin, p['dt_spawn'])
        dteggmax = (p['max_genlength_years'] + 0.5) * 365
        dteggmax = min(dteggmax, forcing['t'][-1])
        dtegg = np.arange(dteggmin, dteggmax + .1, p['dt_spawn'])
    
    ## STRATEGY DICTIONARY
    s = {
        'tdia_exit': tdia_exit,
        'tdia_enter': tdia_enter,
        'dtegg': dtegg
    }
    # Constructing grids using broadcasting
    s['tdia_exit'], s['tdia_enter'], s['dtegg'] = np.meshgrid(tdia_exit, tdia_enter, dtegg, indexing='ij') 
    
    return t0, s