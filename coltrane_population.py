# -*- coding: utf-8 -*-

'''
Coltrane - coltrane_population

@author Lucie Bourreau & Neil Banas
@date 2023/10/26
'''

from tqdm import tqdm
import numpy as np

import warnings
warnings.filterwarnings('ignore')

from timing_combinations import timing_combinations
from add_strategy_to_params import add_strategy_to_params
from coltrane_integrate import coltrane_integrate
from drop_time_series import drop_time_series


def coltrane_population(forcing,p,nargout):
    '''
    v2.0 of the Coltrane model. This has diverged significantly from the
    Coltrane 1.0 model in Banas et al., 2016.
    
    Relationship with the published Coltrane model:
    * the phi model is hereby abandoned
    * R,S have been replaced by R and W = S + R. Allometric formulas using S
        in the paper now use W.
    * the state variables have been largely separated so that not all parts
        require iteration in time. This makes the model cleanly 
        hierarchical, so that one can derive predictions about
            1) development alone (D),
            2) then size and time evolution in surviving cohorts (D,W,R),
            3) then mortality, survivorship, and population dynamics (D,W,R,N,E).
        This will also make it possible for N to be density-dependent in a 
        future version.
    * The myopic criterion for diapause has been replaced by a matrix of entry
        and exit dates, which are analyzed in a brute-force way parallel
        to spawning date.
    * tegg has been replaced by dtegg, which is similar to (tegg - t0).
    
        t           t0          tdia_exit   tdia_enter  dtegg
        timestep    spawn date  exit date   entry date  egg prod date
        (calendar)  (calendar)  (yearday)   (yearday)   (relative to t0)
    
    the last three of these are folded into a single strategy vector s.

    To do:
    Option to run only particular indices into s? useful if saving full time series for
    a small set of optimal cases
    
    Parameters
    ----------
    forcing: dict
        Set of forcing, composed of prey and temperature (surface & deep) cycle over several years.
    p: dict
        Set of parameters (from coltrane_params.py)
    nargout: int
        Number of output arguments you want from the function: 1 for pop, 2 for (pop,popts).
    
    Returns
    -------
    pop: dict
        Contains scalar summaries of what happened in each cohort/strategy combo.
    popts: dict
        Contains full time series of state variables, dF1, and so on.
    '''
    
    # If popts is requested as an output (i.e., nargout=2), we save full time series (which get really big),
    # otherwise discard them and save only summaries
    retain_time_series = True if nargout >= 2 else False

    ts_always_keep = ['t', 'dF1']  # Time series to always save as we go along
    ts_always_omit = list(forcing.keys()) + ['yday']  # Time series to always omit because there are 
                                                      # in the forcing and not generated by the model

    NT = len(forcing['t'])  # Number of timesteps
    t0,s = timing_combinations(forcing, p)
    strategy_fields = list(s.keys())
    NC = len(t0)
    NS = np.prod(s[strategy_fields[0]].shape)
    
    # Run one strategy at a time ---------------------------------------------------------------------------
    
    out = [None] * NS
    print(f"{NT} timesteps x {NC} compupods x {NS} strategies")

    for i in range(NS):
    #for i in tqdm(range(NS), desc="Running strategies"):
        pii = add_strategy_to_params(p, s, i)
        out[i] = coltrane_integrate(forcing, pii, t0) 
        if not retain_time_series:
            out[i] = drop_time_series(out[i], ts_always_keep)
    
    # Clean up output --------------------------------------------------------------------------------------       
    
    # Figure out which strategies produced successful cases

    pop = {}
    pop['level'] = np.zeros((NC,NS))

    for i in range(NS):
        pop['level'][:, i] = out[i]['level']

    f = np.where(np.any(pop['level'] > 0, axis=0)) # Strategies with any complete integrations
    if not np.any(f):
        # If there are not any, return only level (and F1, F2 = 0)
        pop['F1'] = np.zeros(pop['level'].shape)
        pop['F2'] = np.zeros(pop['level'].shape)

        return pop # If any strategy worked, the function stop here and return only pop with zero fitness

    # Rearrange the output from all the individual runs into a single structure (pop)
    # Plus a second one for time series variables if we are keeping them (popts)
    
    popts = {}

    fields = list(out[f[0][0]].keys())
    for k in fields:
        example = out[f[0][0]][k]  # Example for this variable
        if example.ndim == 1:  # If it is a summary field it has only one dimension (i.e. NC)
            pop[k] = np.full((example.shape[0], NS), np.nan) # Create the empty key in pop dict
            for i in range(NS):
                if k in out[i] and out[i][k] is not None:               
                    pop[k][:, i] = out[i][k] # Fill this key in pop dict
        elif (k in ts_always_keep) or (retain_time_series and (k not in ts_always_omit)):
            # if it is a time series (i.e. has 2 dimensions (NT,NC)), and we are retaining time series, 
            # and it is not in the always-omit list -- or if it is dF1 or t, which we always save since they are 
            # necessary to compute the two-generation fitness
            popts[k] = np.full(example.shape + (NS,), np.nan)
            for i in range(NS):
                if k in out[i] and out[i][k] is not None:
                    popts[k][:, :, i] = out[i][k]
    
    # Two-generation fitness -------------------------------------------------------------------------------              
    
    # This calculation is why we need to run cohorts for multiple years (longest lifecycle 
    # plus one) just to evaluate the fitness of the cohorts born in the first year. If you
    # are running a steady-state annual cycle, it would be a big speedup to reimplement a
    # version of the cyclical calculation from v1.0 (as in the 2016 paper).
    
    if 'dF1' in popts:
        F1expected = np.nanmax(pop['F1'], axis=1)  # Expected LEP for each t0, assuming that 
                                                   # the offspring will take the optimal strategy
        F1ex_ = np.interp(popts['t'][:,:, 0], pop['t0'][:, 0], F1expected)
        F1ex_[np.isnan(F1ex_)] = 0
        F1ex_ = np.repeat(F1ex_[:,:,np.newaxis], NS, axis=2)
        dF2 = popts['dF1'] * F1ex_ # Contribution to two-generation fitness at each (t,t0,s)
        pop['F2'] = np.sum(dF2, axis=0)  # Two-generation fitness at each (t0,s)
    
    # Final clean up -----------------------------------------------------------------------------------------         
    
    # Include the strategy fields in pop. (Since coltrane_integrate sees them as model 
    # parameters, they are not saved in the output structures manipulated above).
    
    for k in strategy_fields:
        pop[k] = np.tile(s[k].flatten(), (NC, 1))

    if retain_time_series:
        return pop, popts
    else:
        return pop
