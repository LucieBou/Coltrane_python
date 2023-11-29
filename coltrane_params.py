# -*- coding: utf-8 -*-

'''
Coltrane - coltrane_params

@author Lucie Bourreau & Neil Banas
@date 2023/08/21
'''

def set_default(p0, name, val):
    '''
    Allow user-specified parameters to be inserted into the sequence at the same place they would be declared by
    default; this makes the dependencies among parameters a little more predictable. It may not be necessary at all.
    
    Parameters
    ----------
    p0: dict
        Initial set of parameters (with user-specified parameters empty).
    name: character
        Name of the parameter.
    val: int or char
        Value of the parameter.
        
    Returns
    -------
    p: dict
        Final set of parameters (with the "name" parameters filled).
    '''
    p = p0.copy()
    if name not in p or p[name] is None:
        p[name] = val
    return p

def coltrane_params(**kwargs):
    '''
    Returns a complete set of parameters to run the Coltrane model (coltrane_integrate.ipynb, 
    coltranePopulation.ipynb). Specify which ever non-default values you like and this will fill in the rest.
    
    Parameters
    ----------
    **kwargs: keyworded arguments
    
    Returns
    -------
    p: dictionnary
        Complete set of parameters to run the Coltrane model.
    '''
    
    p = kwargs.copy()
    
    ## LIFE CYCLE
    p = set_default(p, 'dt_spawn', 10)  # resolution (d) of spawning date cases
    p = set_default(p, 'tdia_exit', [])  # set of diapause exit dates to consider
    # if this is empty, constructs a set using dt_dia
    p = set_default(p, 'tdia_enter', [])  # set of diapause entry dates to consider
    # if this is empty, constructs a set using dt_dia
    p = set_default(p, 'dt_dia', 20)
    p = set_default(p, 'dtegg', [])
    # if this is empty, constructs a set using min_ and max_genlength_years
    p = set_default(p, 'min_genlength_years', 0)
    p = set_default(p, 'max_genlength_years', float('inf'))
    # range of generation lengths to evaluate (in integer years)
    
    ## PREDATION
    p = set_default(p, 'preySatVersion', 'default')
    p = set_default(p, 'KsIA', 0.2)
    p = set_default(p, 'iceToSat', 1)
    p = set_default(p, 'Ks', 1.2)  # prey half-saturation; same units as P
    # from tuning, June 2021
    
    ## METABOLISM
    p = set_default(p, 'r_ea', 0.013)  # egg:adult weight ratio if exp_ea=1
    p = set_default(p, 'exp_ea', 0.62)  # egg weight = r_ea * adult weight^exp_ea
    # Kiorboe and Sabatini 1995 (Table 1, Fig 1):
    # sac spawners: r_ea = 0.014, exp_ea = 1
    # broadcast spawners: r_ea = 0.013, exp_ea = 0.62
    p = set_default(p, 'theta', 0.7)  # metabolic scaling exponent
    p = set_default(p, 'u0', 0.008)  # food-saturated development rate at T = 0
    p = set_default(p, 'I0', 0.37)  # food-saturated ingestion rate at S=1 µgC, T=0
    # This was 0.4 in Banas et al. 2016 and Hobbs et al. 2020,
    # but tuning in June 2021 showed that lower I0 was better for the Wa-u0 relationship
    p = set_default(p, 'GGE_nominal', 0.33)  # for relating adult size to I0 and u0.
    p = set_default(p, 'Q10g', 2.5)  # Q10 for growth and ingestion
    p = set_default(p, 'Q10d', 3)  # Q10 for development
    p = set_default(p, 'Df', 0.10)  # age of first feeding (0.10 = start of N3)
    p = set_default(p, 'Ds', 0.35)  # age at which to start storing lipids
    # (0.35 = start of C1)
    p = set_default(p, 'Ddia', 0.6)  # minimum diapause-capable stage (0.6 = start of C4;
    # (alternatively this can be set to 0 or Ds and evaluated in postprocessing)
    p = set_default(p, 'requireActiveSpawning', 0)
    # omit cases in which t0 falls during the diapause period, and set E=0
    # during diapause as well
    p = set_default(p, 'maxReserveFrac', 1)
    # energy gain is discarded if it would push R/W above this level.
    # not an elegant solution, and not a carefully chosen value
    p = set_default(p, 'r_assim', 0.67)  # assimilation efficiency of ingestion
    p = set_default(p, 'rm', 0.8 * 0.17)
    # active metabolism as fraction of max assimilation
    # 0.8*0.17 means that GGE = 0 at P = 0.25 Ks
    p = set_default(p, 'rb', 0.25)  # metabolism at a=0 as fraction of metabolism at a=1
    p = set_default(p, 'rstarv', 0.1)  # starvation tolerance
    
    ## MORTALITY
    # set mortality (m0 = mortality at W = 1 µgC, T = 0)
    p = set_default(p, 'm0_over_GGE_I0', 0.67)
    p = set_default(p, 'm0', p['m0_over_GGE_I0'] * p['GGE_nominal'] * p['I0'])
    p['m0_over_GGE_I0'] = p['m0'] / p['GGE_nominal'] / p['I0']
    # this is probably overcomplicated at this point. Tuning June 2021 says m0 = 0.065

    return p