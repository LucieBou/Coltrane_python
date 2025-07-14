# -*- coding: utf-8 -*-

'''
Coltrane - D_to_stage

@author Lucie Bourreau & Neil Banas
@date 2023/08/24
'''

import numpy as np

def D_to_stage(D):
    '''
    Convert D (Development) to stage, based on schedule in Campbell et al., 2001 (Maps et al., 2014 is similar)
    
    Parameters
    ----------
    D: float
        Development stage between 0 (spawning) and 1 (maturity) of the individual.
    
    Returns
    -------
    stage: int
        Associated stage (between 1 and 13 i.e., E, N1-6, C1-6) of the individual.
    '''
    
    Bele_a = np.array([0, 595, 983, 1564, 2951, 3710, 4426, 5267, 6233, 7370, 8798, 10964, 15047])
    Dstage = Bele_a / Bele_a[-1] # Age of entry into each stage
    
    stage = np.full_like(D, np.nan, dtype=np.float64)
    
    for i in range(1, len(Dstage)):
        stage[(D >= Dstage[i - 1]) & (D < Dstage[i])] = i
    
    stage[D >= 1] = len(Dstage)
    
    return stage