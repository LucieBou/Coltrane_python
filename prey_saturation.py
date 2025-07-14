# -*- coding: utf-8 -*-

'''
Coltrane - prey_saturation

@author Lucie Bourreau & Neil Banas
@date 2023/08/23
'''

import numpy as np

def prey_saturation(v0, p):
    '''
    Calculate prey saturation based on some settings within p, which is the structure that comes out of 
    coltrane_params.py. In Coltrane 1.0, all the complexities of this were handled in coltraneForcing.m, 
    before the model was run, but doing things in this order makes it possible to vary assumptions about the 
    forcing as part of a big parameterisation experiment.
    This function is adapted to different forcing. For example, prey cycles can be composed of diatoms, 
    flagellates, etc.
    
    Parameters
    ----------
    v0: dict
        Set of forcing. Like the forcing input in timing_combinations.py.
    p: dict
        Set of parameters (from coltrane_params.py).
    
    Returns
    -------
    v: dict
        Set of forcing with the additionnal prey saturation parameter 'sat'.
    '''
    
    v = v0.copy()

    if p['preySatVersion'].lower() == 'biomas_dia21':
        if 'Ptot' not in v:
        # Combined the different prey categories into one
            v['Ptot'] = v['flagel'] + v['diatom']
            
        # Prey saturation considering water-column prey only
        # To be able to do arithmetic operations on dict we need to convert the list from the dict to numpy array
        v['satWC'] = np.array(v['Ptot']) / (p['Ks'] + np.array(v['Ptot']))
        
        # Ice alage index (IAind) mimicking Castellani et al., 2017, a function of yearday and latitude multiplied 
        # by ice cover        
        t_init = np.maximum(45, 2.78 * np.array(v['y']) - 117) # v['y'] corresponds to latitude
        t_max = np.maximum(45, 2.08 * np.array(v['y']) - 30) # linear fits to Castellani et al., 2017 Table 5
        t_end = 200 # mid July
        IAind = np.zeros_like(v['y'])
        f = np.where((v['yday'] >= t_init) & (v['yday'] <= t_max))[0]
        IAind[f] = np.array(v['ice'])[f] * (np.array(v['yday'])[f] - t_init[f]) / (t_max[f] - t_init[f])
        f = np.where((v['yday'] >= t_max) & (v['yday'] <= t_end))[0]
        IAind[f] = np.array(v['ice'])[f] * (1 - (np.array(v['yday'])[f] - t_max[f]) / (t_end - t_max[f]))
        
        # Prey saturation considering ice-algae only, based on an effective
        # half-saturation (in index units, not chl units). I briefly tried an adjustable
        # iceToSat multiplying this expression but tuning suggested that 1 is about right-- whereas KsIA << 1
        v['satIA'] = IAind / (p['KsIA'] + IAind)
        
        v['sat'] = np.maximum(v['satWC'], v['satIA'])
        
    elif p['preySatVersion'].lower() == 'biomas_dia19a':
                # Same as biomas_dia19 below, but tIA is a function of latitude chosen so that the time point 1/3
                # of the way through the dtIA interval aligns with the bloom max a function of latitude from 
                # Castellani et al. 2017 so dtIA and iceToSat are the free parameters.

        if 'Ptot' not in v:
            v['Ptot'] = v['flagel'] + v['diatom']

            v['satWC'] = np.array(v['Ptot']) / (p['Ks'] + np.array(v['Ptot']))

            tIA = np.maximum(45, 2.08 * np.array(v['y']) - 30) - p['dtIA'] / 3
            v['satIA'] = (np.array(v['ice']) > 0.15) * (np.array(v['yday']) > tIA) * (np.array(v['yday']) < 200) * p['iceToSat']

            yr = np.ceil((v['t'] - v['t'][0]) / 365)
            for n in range(1, int(np.max(yr)) + 1):
                satn = np.zeros_like(v['satIA'])
                satn[yr == n] = v['satIA'][yr == n]
                satn[np.cumsum(satn.astype(int)) > p['dtIA'] / (v['t'][1] - v['t'][0])] = 0
                v['satIA'][yr == n] = satn[yr == n]

                v['sat'] = np.maximum(v['satWC'], v['satIA'])

    elif p['preySatVersion'].lower() == 'biomas_dia19':
        if 'Ptot' not in v:
            v['Ptot'] = v['flagel'] + v['diatom']

            v['satWC'] = np.array(v['Ptot']) / (p['Ks'] + np.array(v['Ptot']))

            v['satIA'] = (np.array(v['ice']) > 0.15) * (np.array(v['yday']) > p['tIA']) * (np.array(v['yday']) < 200) * p['iceToSat']

            yr = np.ceil((v['t'] - v['t'][0]) / 365)
            for n in range(1, int(np.max(yr)) + 1):
                satn = np.zeros_like(v['satIA'])
                satn[yr == n] = v['satIA'][yr == n]
                satn[np.cumsum(satn.astype(int)) > p['dtIA'] / (v['t'][1] - v['t'][0])] = 0
                v['satIA'][yr == n] = satn[yr == n]

                v['sat'] = np.maximum(v['satWC'], v['satIA'])

    elif p['preySatVersion'].lower() == 'satellite_dia18':
        # Prey saturation considering water-column prey only 
        v['chl'][np.isnan(v['chl']) & (np.array(v['ice']) > 0.15)] = p['chlUnderIce']
        v['chl'][np.isnan(v['chl'])] = p['chlUnderPersistentCloud']
        v['satWC'] = np.array(v['chl']) / (p['Ks'] + np.array(v['chl']))

        # Prey saturation considering (a guess at) ice algae only. Days that ice cover is > 15%, 
        # after yearday _tIA_ and before yearday (365-tIA), all weighted by a highly uncertain 
        # weighting factor iceToSat
        v['satIA'] = (np.array(v['ice']) > 0.15) * (np.array(v['yday']) > p['tIA']) * (np.array(v['yday']) < 270) * p['iceToSat']

        v['sat'] = np.maximum(v['satWC'], v['satIA'])

    else:
        v['sat'] = np.array(v['P']) / (p['Ks'] + np.array(v['P']))
    
    return(v)