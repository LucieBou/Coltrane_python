#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Coltrane - Forcing time series

@author: luciebourreau
@date 2023/11/29
"""

import os
import numpy as np
import pandas as pd
from Create_ts_BB_or_NOW_v2 import time_series_NOW_BB

def coltrane_forcing(region, Nyears):
    """
    Create the forcing time series necessary to run the Coltrane model 
    from a given Arctic region.
    
    Disko Bay: one seasonal cycle of temperature and 
    phytoplankton and microzooplankton prey, based on the 1996–1997 time series 
    described by Madsen et al. (2001). We use this particular dataset not 
    primarily as a guide to the current or future state of Disko Bay but rather 
    as a specific circumstance in which the life-history patterns of three 
    coexisting Calanus spp. (C. finmarchicus, C. glacialis, C. hyperboreus) 
    were documented (Madsen et al., 2001).
    
    TO BE UPDATED WITH OTHER REGIONS FORCING.


    Parameters
    ----------
    region: string
        Region where the forcing are coming from.
    Nyears: scalar
        Number of years to run the model.

    Returns
    -------
    forcing: dict
        Forcing time series.
    """
    
    forcing = {}
    
    if region == "DiskoBay": 
        
        forcing['t'] = np.arange(0, 365)  # start with one year
        forcing['T0'] = np.nan * forcing['t']
        forcing['Td'] = np.nan * forcing['t']
        forcing['P'] = np.nan * forcing['t']

        ## Surface temperature

        Tjan1 = -0.7  # T on Jan 1
        Tmin = -1.8  # winter min T
        Tmax = 3.7  # summer max T
        tTmin = 105  # time when spring T increase starts
        tTmax = 250  # time of T max
        nn = np.where((forcing['t']+1 >= 0) & (forcing['t']+1 <= tTmin))[0]
        forcing['T0'][nn] = np.linspace(Tjan1, Tmin, len(nn))
        nn = np.where((forcing['t']+1 >= tTmin) & (forcing['t']+1 <= tTmax))[0]
        forcing['T0'][nn] = np.linspace(Tmin, Tmax, len(nn))
        nn = np.where((forcing['t']+1 >= tTmax) & (forcing['t']+1 <= 365))[0]
        forcing['T0'][nn] = np.linspace(Tmax, Tjan1, len(nn))

        ## Deep temperature

        Td0 = 1  # deep T year-round
        forcing['Td'] = np.full_like(forcing['t'], Td0)

        ## Chlorophyll (prey)

        P0win = 0.1
        P0spr = 13
        P0sum = 0.05 * P0spr
        P0aut = 5
        dtPspr = 15  # bloom duration in spring
        dtPaut = 30 # bloom duration in autumn
        tPspr = 150 # time bloom initiation in spring
        tPaut = 225 # time bloom initiation in autumn
        forcing['P'] = P0win * np.ones_like(forcing['t'])
        forcing['P'] = np.maximum(forcing['P'], P0spr * np.exp(-((forcing['t']+1 - tPspr) / dtPspr) ** 2))
        forcing['P'] = np.maximum(forcing['P'], P0aut * np.exp(-((forcing['t']+1 - tPaut) / dtPaut) ** 2))
        fsum = (forcing['t']+1 > tPspr) & (forcing['t']+1 < tPaut)
        forcing['P'][fsum] = np.maximum(forcing['P'][fsum], P0sum)

        ## Repeat this cycle for Nyears

        forcing['t'] = np.arange(0, Nyears * 365)
        forcing['T0'] = np.tile(forcing['T0'], Nyears)
        forcing['Td'] = np.tile(forcing['Td'], Nyears)
        forcing['P'] = np.tile(forcing['P'], Nyears)
        
    if region == "NOW":
        
        forcing = time_series_NOW_BB("NOW",
                                     2013,
                                     2013,
                                     #"/Users/luciebourreau/Library/CloudStorage/OneDrive-UniversitéLaval/PhD_ULaval/Data/NOW_BB_profiles_Inge/",
                                     True)
        
        forcing['t'] = np.arange(0, Nyears * 365)
        
        keys = list(forcing.keys())
        keys.remove('t')
        keys.remove('year')
        
        for key in keys:
            forcing[key] = np.tile(forcing[key], Nyears)
        
    if region == "BB":
        
        forcing = time_series_NOW_BB("BB",
                                     2013,
                                     2013,
                                     #"/Users/luciebourreau/Library/CloudStorage/OneDrive-UniversitéLaval/PhD_ULaval/Data/NOW_BB_profiles_Inge/",
                                     True)
        
        forcing['t'] = np.arange(0, Nyears * 365)
        
        keys = list(forcing.keys())
        keys.remove('t')
        keys.remove('year')
        
        for key in keys:
            forcing[key] = np.tile(forcing[key], Nyears)
    
    if region == "Qik_mod_2015":

        forcing['t'] = np.arange(0, 365)  # start with one year
        forcing['T0'] = np.nan * forcing['t']
        forcing['Td'] = np.nan * forcing['t']
        forcing['P'] = np.nan * forcing['t']

        # Load the data
        script_dir = os.path.dirname(os.path.abspath(__name__))  # Répertoire du script
        NEMO_file_path = os.path.join(script_dir, 'model/NEMO_T0_Td_Qik_2016.csv')
        fluo_file_path = os.path.join(script_dir, 'model/max_fluo_mod_obs_qik_2015_lowess016.csv')

        NEMO = pd.read_csv(NEMO_file_path)
        max_fluo = pd.read_csv(fluo_file_path)

        forcing['T0'] = NEMO['NEMO_T0']
        forcing['Td'] = NEMO['NEMO_Td']
        forcing['P'] = max_fluo['Max_fluo_obs_mod']

        ## Repeat this cycle for Nyears

        forcing['t'] = np.arange(0, Nyears * 365)
        forcing['T0'] = np.tile(forcing['T0'], Nyears)
        forcing['Td'] = np.tile(forcing['Td'], Nyears)
        forcing['P'] = np.tile(forcing['P'], Nyears)

    if region == "Qik_obs_2015":

        forcing['t'] = np.arange(0, 365)  # start with one year
        forcing['T0'] = np.nan * forcing['t']
        forcing['Td'] = np.nan * forcing['t']
        forcing['P'] = np.nan * forcing['t']

        # Load the data
        # script_dir = os.path.dirname(os.path.abspath(__name__))  # Répertoire du script
        # NEMO_file_path = os.path.join(script_dir, 'model/NEMO_T0_Td_Qik_2016.csv')
        # fluo_file_path = os.path.join(script_dir, 'model/max_fluo_obs_qik_2015_lowess025.csv')

        script_dir = os.path.dirname(os.path.abspath(__file__))  # Script repository
        NEMO_file_path = os.path.join(script_dir, 'NEMO_T0_Td_Qik_2016.csv')
        fluo_file_path = os.path.join(script_dir, 'max_fluo_obs_qik_2015_lowess025.csv')


        NEMO = pd.read_csv(NEMO_file_path)
        max_fluo = pd.read_csv(fluo_file_path)

        forcing['T0'] = NEMO['NEMO_T0']
        forcing['Td'] = NEMO['NEMO_Td']
        forcing['P'] = max_fluo['Max_fluo_obs_fill']

        ## Repeat this cycle for Nyears

        forcing['t'] = np.arange(0, Nyears * 365)
        forcing['T0'] = np.tile(forcing['T0'], Nyears)
        forcing['Td'] = np.tile(forcing['Td'], Nyears)
        forcing['P'] = np.tile(forcing['P'], Nyears)

    return forcing
