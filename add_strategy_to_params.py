# -*- coding: utf-8 -*-

'''
Coltrane - add_strategy_to_params

@author Lucie Bourreau & Neil Banas
@date 2023/10/02
'''

def add_strategy_to_params(p,s,ind):
    '''
    Copy a single value of the strategy vector s into p as parameters.
    
     Parameters
    ----------
    p: dict
        Set of parameters.
    s: dict
        Timing strategy vector: diapause exit date, diapause entry date and egg production range.
    ind: float
        Index, i.e. strategy number.
    
    Returns
    -------
    pii: dict
        Set of parameters with only one strategy.
    '''
    
    pii = p.copy()
    fields = s.keys()
    for k in fields:
        sflat = s[k].flatten()
        pii[k] = sflat[ind]
        
    return pii