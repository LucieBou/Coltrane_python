# -*- coding: utf-8 -*-

'''
Coltrane - coltrane_integrate

@author Lucie Bourreau & Neil Banas
@date 2023/09/23
'''

# Load the packages

import numpy as np
from datetime import date 

# Load the functions

from yearday import yearday
from prey_saturation import prey_saturation
from stage_to_D import stage_to_D

def coltrane_integrate(forcing,p,t0):
    '''
    Calculates a(t), D(t), W(t), R(t), N(t), E(t), and dF1 for a set of spawning dates t0 and a fixed timing
    strategy (s), along with a variety of summary metrics.
    The three fields of s (tdia_exit, tdia_enter, and dtegg) should be inside p as additional scalar parameters.
    
    This is where most of the actual model equations live. coltrane_integrate performs the actual time 
    integration for a range of t0 values but a single timing strategy, whose fields (tdia_enter, etc) are 
    copied into the parameter list p. The fields of v are a mix of full time series of size [NT NC] 
    and summary metrics of size [1 NC].
    
    Parameters
    ----------
    forcing: dict
        Set of forcing, composed of prey and temperature (surface & deep) cycle over several years.
    p: dict
        Set of parameters (from coltraneParams.ipynb) with the addition of a single value of the strategy vector s 
        (from addStrategyToParams function in coltrane_community.py): the three fields of s (tdia_exit, tdia_enter, 
        and dtegg) should be inside p as additional scalar parameters.
    t0: array
        Spawning dates (from timing_combinations.py).
    
    Returns
    -------
    v: dict
        Time series of state variables and summary metrics for one timing strategy.
        
    '''
    
    # Prepare the outputs -----------------------------------------------------------------------------------       
    
    v = forcing.copy()
    v['t0'] = t0.flatten()
    NC = v['t0'].shape[0]  # number of cohorts (t0 values)
    NT = v['t'].shape[0]   # number of timesteps
    

    # Extend the forcing NC times in the v dict
    for key in forcing:
        v[key] = np.transpose(np.tile(v[key], (NC,1)))
    
    # Timebase
    v['yday'] = yearday(v['t'])  # make sure this is filled in and consistent
    dt = v['t'][2,0] - v['t'][1,0]  # determine simulation timestep from the forcing time series 
                                    # (v['t'][0,0] corresponds to the first value of the first array)
                                    # (v['t'][0,1] corresponds to the second value of the first array)

    # Define output variables that exist for all cases, even if there are no valid solutions 
    # to the model integration
    v['F1'] = np.zeros(NC)
    v['level'] = np.zeros(NC)
        # level 0 = logical inconsistency (between development and dtegg)
        # level 1 = successful diapause
        # level 2 = reaches adulthood
        # and past there, fitness provides classifications
    
    # Activity: a(t) ----------------------------------------------------------------------------------------                      
    
    v['a'] = np.double(~((v['yday'] >= p['tdia_enter']) | (v['yday'] <= p['tdia_exit']))) # There is an offset 
    # of 1 with respect to the matlab code and a value of 1 is missing each time. 
    # Is this due to the difference in days?
    isalive = v['t'] >= np.tile(v['t0'], (NT,1))
    
    if p['requireActiveSpawning']:
        # Eliminate t0 values falling during diapause
        t0_yday = yearday(v['t0'])
        activeSpawning = ~((t0_yday >= p['tdia_enter']) | (t0_yday <= p['tdia_exit']))
        isalive = isalive & np.tile(activeSpawning, (NT, 1))

        # Likewise, t0 values in which t0+dtegg falls during diapause (this eliminates some redundancy)
        ta_yday = yearday(v['t0'] + p['dtegg'])
        activeNextSpawning = ~((ta_yday >= p['tdia_enter']) | (ta_yday <= p['tdia_exit']))
        isalive = isalive & np.tile(activeNextSpawning, (NT, 1))
    
    # Temperature response factors: qd, qg  -----------------------------------------------------------------       
    
    v['temp'] = v['a'] * v['T0'] + (1 - v['a']) * v['Td']
    
    qd = np.zeros((NT, NC))
    qd[isalive] = (p['Q10d'] ** 0.1) ** v['temp'][isalive]
    
    qg = np.zeros((NT, NC))
    qg[isalive] = (p['Q10g'] ** 0.1) ** v['temp'][isalive]
    
    # Prey saturation ---------------------------------------------------------------------------------------
    
    v = prey_saturation(v,p)
    if v['sat'].shape[1] == 1:
        v['sat'] = np.tile(v['sat'], (1,NC))
    
    # Development: D(t) -------------------------------------------------------------------------------------              
    
    dDdt = isalive * p['u0'] * qd  # nonfeeding formula
    
    v['D'] = np.cumsum(dDdt, axis=0) * dt
    isfeeding = v['D'] >= p['Df']
    
    dDdt_feeding = isalive * p['u0'] * qd * v['a'] * v['sat']
    dDdt[isfeeding] = dDdt_feeding[isfeeding]  # full formula
    
    v['D'] = np.cumsum(dDdt, axis=0) * dt
    v['D'][v['D'] > 1] = 1
    isfeeding = v['D'] >= p['Df']   # recalculated because cumsum() is one timestep off from a 
                                    # true forward integration
    
    # The life history is divided into two phases, growth and egg production. 
    # This is equivalent to juvenile and adult if the animals delay their
    # final maturation until just before egg prod. Begins, but not if they enter
    # C6 and then delay egg prod. (e.g. in the case of overwintering C6's)
    
    isineggprod = v['t'] >= v['t0'] + p['dtegg']
        
    # Check D(t) for consistency with dtegg
    toolate = np.any((v['D'] < 1) & (isineggprod), axis=0)
    v['D'][:, toolate] = np.nan
    # if np.all(toolate):
    #     return None # if there are no good solutions at all, don't bother with the rest of the model
    
    # Calculate D in middle of first winter, just as a diagnostic
    year = v['t0'] // 365
    first31dec = 365 + year * 365
    first31dec = np.tile(first31dec, (NT,1))

    is31dec = np.abs(v['t'] - first31dec) == np.tile(np.min(np.abs(v['t'] - first31dec), axis=0), (NT,1))
    is31dec = is31dec & np.cumsum(is31dec, axis=0) == 1

    v['D_winter'] = np.reshape(v['D'][is31dec], NC)
    
    # Flag time points at which the animal is in diapause but at a
    # diapause-incapable stage, and mark these cases as dead
    isactive = isalive & ((v['a'] == 1) | (~isfeeding))
    hasbeenactive = np.cumsum(isactive, axis=0) >= 1
    isfailingtodiapause = isalive & hasbeenactive & (v['a'] == 0) & (v['D'] < p['Ddia'])
    hasfailedtodiapause = np.cumsum(isfailingtodiapause, axis=0) > 1
    v['D'][hasfailedtodiapause] = np.nan
    v['level'][~np.any(hasfailedtodiapause, axis=0)] = 1  # level 1 = successful diapause
    
    isalive = isalive & ~np.isnan(v['D'])
    isfeeding = isfeeding & isalive
    isineggprod = isineggprod & isalive

    # Date on which D=1 is reached (another diagnostic)
    tD1 = v['t'].astype(float).copy()
    hasreachedadulthood = (v['D'] < 1) | (~np.isfinite(v['D']))
    tD1[hasreachedadulthood] = np.nan
    v['tD1'] = np.nanmin(tD1, axis = 0)
    v['level'][np.any(v['D'] == 1, axis = 0)] = 2  # level 2 = reaches adulthood

    # Energy gain, Growth and Egg production: G(t), W(t), R(t), E(t) ----------------------------------------       
    
    # Estimation of adult size based on mean forcing temperature (Wa_theo), 
    # and corresponding estimation of egg size based on this temperature (We_theo).
    T_nominal = np.mean(v['temp'], axis=0)
    qoverq = (p['Q10g'] / p['Q10d']) ** (T_nominal / 10)
    co = (1 - p['theta']) * p['GGE_nominal'] * (1 - p['Df']) * qoverq
    v['Wa_theo'] = (co * p['I0'] / p['u0']) ** (1 / (1 - p['theta']))
    v['We_theo'] = p['r_ea'] * v['Wa_theo'] ** p['exp_ea']
    
    # Compute the state variables
    
    v['G'] = np.zeros((NT, NC)) # Growth
    v['W'] = np.tile(v['We_theo'], (NT, 1)) # Weight
    v['R'] = np.zeros((NT, NC)) # Reserves
    v['Einc'] = np.zeros((NT, NC)) # Income breeding
    v['E'] = np.zeros((NT, NC)) # Egg production
    astar = p['rb'] + (1 - p['rb']) * v['a']
    
    for n in range(0,NT-1):
        # For each timestep, we calculate the growth and energy gain for each cohort.
        f = isfeeding[n, :]
        e = isineggprod[n, :]

        # Net gain
        Imax_nf = qg[n, f] * p['I0'] * v['W'][n, f] ** (p['theta'] - 1)
        I_nf = v['a'][n, f] * p['r_assim'] * v['sat'][n, f] * Imax_nf
        M_nf = p['rm'] * astar[n, f] * Imax_nf
        v['G'][n, f] = I_nf - M_nf
        GWdt = v['G'][n, :] * v['W'][n, :] * dt

        # Allocation to growth
        dW = GWdt.copy()
        dW[(dW > 0) & e] = 0  # no storage or growth once egg production has begun
        v['W'][n+1, :] = np.maximum(0, np.nan_to_num(v['W'][n, :] + dW, nan=0))

        # Allocation to reserves
        fr = (v['D'][n, :] - p['Ds']) / (1 - p['Ds']) # Ds: age at which to start storing lipids (around C1)
        fr = np.maximum(0, np.minimum(1, np.nan_to_num(fr, nan=1)))
        fr[GWdt < 0] = 1  # all net losses come from R
        v['R'][n+1, :] = v['R'][n, :] + fr * dW 
        # if R/W > maxReserveFrac, we added too much: throw some of dW away
        excess = np.maximum(0, np.nan_to_num(v['R'][n+1, :] - p['maxReserveFrac'] * v['W'][n+1, :], nan=0))
        v['W'][n+1, :] = v['W'][n+1, :] - excess
        v['R'][n+1, :] = v['R'][n+1, :] - excess

        # Income egg production
        v['Einc'][n, :] = np.maximum(0, np.nan_to_num(GWdt, nan=0)) / dt * e

        # Capital egg production
        Emax = np.zeros_like(GWdt)
        Emax[f] = p['r_assim'] * Imax_nf * e[f] * v['W'][n, f]
        Ecap = np.maximum(0, np.nan_to_num(Emax - v['Einc'][n, :], nan=0))
        dR = np.minimum(np.maximum(0, np.nan_to_num(v['R'][n+1, :], nan=0)), Ecap * dt)
        v['E'][n, :] = v['Einc'][n, :] + dR / dt
        v['W'][n+1, :] = v['W'][n+1, :] - dR
        v['R'][n+1, :] = v['R'][n+1, :] - dR
    
    # Adult size Wa, Ra (= size at the moment egg prod begins)
    last = ~isineggprod[0:-1, :] & isineggprod[1:, :]
    v['Wa'] = np.nanmax(v['W'][0:-1, :] * last, axis=0)
    v['Ra'] = np.nanmax(v['R'][0:-1, :] * last, axis=0)

    # Update the estimate of We
    v['We'] = p['r_ea'] * v['Wa'] ** p['exp_ea']

    # W and R at all stages
    stages = ['N6', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6']
    for i in range(1, len(stages) - 1):
        Dstart_i = 0.5 * (stage_to_D(stages[i - 1]) + stage_to_D(stages[i]))
        Dend_i = 0.5 * (stage_to_D(stages[i]) + stage_to_D(stages[i + 1]))
        atstage_i = np.double(np.logical_and(v['D'] >= Dstart_i, v['D'] < Dend_i) & isalive)
        atstage_i[atstage_i == 0] = np.nan
        v['W_' + stages[i]] = np.nanmean(v['W'] * atstage_i, axis=0)
        v['R_' + stages[i]] = np.nanmean(v['R'] * atstage_i, axis=0)
    
    # W, R for winter late stages
    # (this is not a general definition of winter, but calculating it here avoids the
    # need to save full time-series output)
    D_C5_start = 0.5 * (stage_to_D('C4') + stage_to_D('C5'))
    iswin = np.logical_or(v['yday'] >= date(2000,11,1).timetuple().tm_yday, v['yday'] <= date(2000,2,28).timetuple().tm_yday)
    iswinlatestage = np.double(np.logical_and(np.logical_and(v['D'] >= D_C5_start, isalive), iswin))
    iswinlatestage[iswinlatestage == 0] = np.nan
    v['W_C56win'] = np.nanmean(v['W'] * iswinlatestage, axis=0)
    v['R_C56win'] = np.nanmean(v['R'] * iswinlatestage, axis=0)
    
    # Check for starvation
    isstarving = v['R'] < -p['rstarv'] * v['W']
    isalive = isalive & (np.cumsum(isstarving, axis=0) == 0)
    isineggprod = isineggprod & isalive
    
    v['E'][~isalive] = 0
    v['Einc'][~isalive] = 0
    
    # Capital fraction of egg production
    v['capfrac'] = np.nansum(v['E'] - v['Einc'], axis=0) / np.nansum(v['E'], axis=0)
    
    # Mortality and survivorship: N(t) ---------------------------------------------------------------------
    
    v['m'] = np.zeros((NT, NC))
    v['m'][isalive] = p['m0'] * qg[isalive] * v['a'][isalive] * v['W'][isalive]**(p['theta'] - 1) # mort. rate at T, size
    v['lnN'] = np.cumsum(-v['m'], axis=0) * dt
    # Calculate adult recruitment (= recruitment at the moment egg prod begins)
    lnNa = v['lnN'].copy()
    lnNa[~isineggprod] = np.nan
    v['Na'] = np.exp(np.nanmax(lnNa, axis=0))
    
    # Contributions to fitness at each t -------------------------------------------------------------------       
    
    v['dF1'] = np.real(v['E'] * np.exp(v['lnN']) * dt) / v['We_theo']
    v['dF1'][np.isnan(v['dF1'])] = 0
    v['F1'] = np.sum(v['dF1'], axis=0)
    
    # Timing metrics ---------------------------------------------------------------------------------------       
    
    v['tEcen'] = np.sum(v['t'] * v['dF1'], axis=0) / v['F1']
    # tEcen - t0 is generation length: similar to, but more accurate than, dtegg - t0
    
    GWNdt = v['G'] * v['W'] * np.exp(v['lnN']) * dt
    GWNdt[np.isnan(GWNdt)] = 0
    GWNdt = np.maximum(0, GWNdt)
    v['tGain'] = np.sum(v['t'] * GWNdt, axis=0) / np.sum(GWNdt, axis=0)  # Center of mass of gain G*W*N

    mWNdt = v['m'] * v['W'] * np.exp(v['lnN']) * dt
    mWNdt[np.isnan(mWNdt)] = 0
    v['tYield'] = np.sum(v['t'] * mWNdt, axis=0) / np.sum(mWNdt, axis=0)  # Center of mass of yield m*W*N

    mRNdt = v['m'] * v['R'] * np.exp(v['lnN']) * dt
    mRNdt[np.isnan(mRNdt)] = 0
    v['tYieldR'] = np.sum(v['t'] * mRNdt, axis=0) / np.sum(mRNdt, axis=0)  # Center of mass of lipid yield m*R*N

    mWNdt = v['m'] * v['W'] * np.exp(v['lnN']) * dt * (v['D'] >= D_C5_start)
    mWNdt[np.isnan(mWNdt)] = 0
    v['tYieldC56'] = np.sum(v['t'] * mWNdt, axis=0) / np.sum(mWNdt, axis=0)  # Center of mass of C5-6 yield

    # Scalar metrics from forcing --------------------------------------------------------------------------       
    
    fields = ['Ptot', 'T0', 'Td', 'sat', 'ice', 'satIA', 'satWC']
    a0 = v['a'].copy()
    a0[~np.isfinite(a0)] = 0

    for k in range(len(fields)):
        if fields[k] in v:
            f0 = v[fields[k]].copy()
            f0[~np.isfinite(f0)] = 0
            if f0.shape[1] == 1:
                f0 = np.tile(f0, (1, NC))
            v[fields[k] + '_avg'] = np.nanmean(f0, axis=0)
            v[fields[k] + '_active'] = np.sum(f0 * a0, axis=0) / np.sum(a0, axis=0)
            isgrowing = isalive & ~isineggprod
            v[fields[k] + '_growing'] = np.sum(f0 * isgrowing, axis=0) / np.sum(isgrowing, axis=0)

    # Clean up ---------------------------------------------------------------------------------------------        
    # Blank out nonliving portions of the time series
    
    v['a'][~isalive] = np.nan
    v['D'][~isalive] = np.nan
    v['W'][~isalive] = np.nan
    v['R'][~isalive] = np.nan
    v['lnN'][~isalive] = -np.inf
    v['E'][~isalive] = 0
    
    return v
