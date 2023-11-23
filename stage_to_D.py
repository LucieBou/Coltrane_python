# -*- coding: utf-8 -*-

'''
Coltrane - stage_to_D

@author Lucie Bourreau & Neil Banas
@date 2023/08/23
'''

def stage_to_D(stage):
    '''
    Quick lookup for middle of each developmental stage as a fraction of total development, based on schedule 
    in Campbell et al. 2001.
    
    Parameters
    ----------
    stage: chr
        Stage of the individual (N1-6 or C1-6).
    
    Returns
    -------
    D: float
        Associated development stage (between 0 (spawning) and 1 (maturity)).
    '''
    
    if stage.upper() == 'C6':
        return 1
    else:
        stages = ['E', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'C1', 'C2', 'C3', 'C4', 'C5']
        Bele_a = [0, 595, 983, 1564, 2951, 3710, 4426, 5267, 6233, 7370, 8798, 10964, 15047]
        n = stages.index(stage.upper())
        return 0.5 * (Bele_a[n + 1] + Bele_a[n]) / Bele_a[-1]