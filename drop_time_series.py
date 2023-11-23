# -*- coding: utf-8 -*-

'''
Coltrane - drop_time_series

@author Lucie Bourreau & Neil Banas
@date 2023/10/31
'''

def drop_time_series(out, ts_alwaysKeep):
    '''
    Remove time series variables to save (a massive amount of) memory
    
     Parameters
    ----------
    out: dict
        Compilation of the computed state variables in coltrane_integrate.py for a given strategy.
    ts_alwaysKeep: list
        Time series to keep.
    
    Returns
    -------
    out1: dict
        Same as out but without the time series except the ones to keep.
        
    '''
    out1 = {}
    fields = out.keys()

    for k in fields:
        if out[k].ndim == 1:
            out1[k] = out[k]

    for k in ts_alwaysKeep:
        if k in out:
            out1[k] = out[k]
            
    return out1
