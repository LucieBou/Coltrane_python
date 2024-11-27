#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function to create time series from NOW or BB profiles

@author: Lucie Bourreau
@date: 2024/03/29
"""

import os
import pickle
import numpy as np

def time_series_NOW_BB(region, year1, year2, extend):
    """
    Time series of different environmental and biological variables from a 
    coupled biogeochemical model on an Arctic and Northern Hemisphere Atlantic 
    1/4° resolution (ANHA4) tripolar grid configuration.
    
    BB and NOW were delimited thanks to a mask and the profiles per day correspond
    to the mean values of the profiles in this area.
    
    The time series are composed of:
        - 't_yr': time in days (from (0 to 365)* nb_years)
        - 't': time in days (from 0 to nb_years * 365)
        - 'year': years of the time serie
        - 'T0': Surface température (in ° Celsius) - mean value between the surface and 100m depth
        - 'Td': Deep temperature (in ° Celsius) - mean value below 300m depth
        - 'P': Prey (in mg.chl.m-3) - maximum of chlorophylle in the profile (from the sum of diatoms and flagellates)
        - 'diat': Diatoms (in mg.chl.m-3) - maximum of chlorophylle in the profile
        - 'flag': Flagellates (in mg.chl.m-3) - maximum of chlorophylle in the profile
        - 'microz': microzooplankton (in mmolN.m-3) - mean values within the upper 100m depth
        - 'mesoz': mesozooplankton (in mmolN.m-3) - mean values within the upper 100m depth
        - 'ialg': ice algae (in mmolN.m-3) - maximum in the profile

    Parameters
    ----------
    region : str
        'NOW' or 'BB' i.e., North Water Polynya or Baffin Bay.
    year1 : int
        Start year of time series.
    year2 : int
        End year of time series.
    folder_path : str
        Path of the folder where the model's outputs are store.
    extend : boolean
        True or False 
            - True if you want to have a value per day (linear interpolation)
            - False if you want the original values, i.e., every 5 days.

    Returns
    -------
    ts : dict
        Time series.

    """
    
    #os.chdir(folder_path)
    
    # Initialize the ts outputs
    ts = {'t_yr': [],
          't': [],
          'year': [],
          'T0': [],
          'Td': [],
          'P': [],
          'diat': [],
          'flag': [],
          'microz': [],
          'mesoz': [],
          'ialg': []}
        
    if region == 'BB':
        
        # Open the dict with all the profiles (all years)

        script_dir = os.path.dirname(os.path.abspath(__file__))  # Répertoire du script
        file_path = os.path.join(script_dir, 'BB_profiles_2002-2021_dict.pkl')
        
        with open(file_path, 'rb') as file:
                BB_profiles_all = pickle.load(file)
        
        # Select only the years we want
        if year1 == year2:
            BB_profiles = BB_profiles_all[year1]
            
            # Time
            ts['t'] = list(range(0,365,5))
            ts['t_yr'] = list(range(0,365,5))
            ts['year'] = list(range(year1, year2+1))
            
            # Surface Temperature (° Celsius)
            temp = BB_profiles['TEMsurf']
            depth = BB_profiles['depth']
            
            ts['T0'] = np.nanmean(np.where(depth <= 100, temp, np.nan), axis=0)
            
            # Deep Temperature (° Celsius)
            ts['Td'] = np.nanmean(np.where(depth >= 300, temp, np.nan), axis=0)
            
            # Prey (mg.chl.m-3)
            prey = BB_profiles['dia2prof'] + BB_profiles['fla2prof']
            
            ts['P'] = np.nanmax(prey, axis=0)
            
            # Diatoms (mg.chl.m-3)
            diat = BB_profiles['dia2prof']
            ts['diat'] = np.nanmax(diat, axis=0)
            
            # Flagellates (mg.chl.m-3)
            flag = BB_profiles['fla2prof']
            ts['flag'] = np.nanmax(flag, axis=0)
            
            # Microzooplankton (mmolN.m-3)
            microz = BB_profiles['micprof']
            ts['microz'] = np.nanmean(np.where(depth <= 100, microz, np.nan), axis=0)
            
            # Mesozooplankton (mmolN.m-3)
            mesoz = BB_profiles['mesprof']
            ts['mesoz'] = np.nanmean(np.where(depth <= 100, mesoz, np.nan), axis=0)
            
            # Ice Algae (mmolN.m-3)
            ialg = BB_profiles['ialgprof']
            ts['ialg'] = np.nanmax(ialg, axis=0)
            
            if extend:
                
                ts = {'t_yr': list(range(0,365)),
                      't': list(range(0,365)),
                      'year': list(range(year1, year2+1)),
                      'T0': np.interp(np.arange(365), np.arange(0, 365, 5), ts['T0']),
                      'Td': np.interp(np.arange(365), np.arange(0, 365, 5), ts['Td']),
                      'P': np.interp(np.arange(365), np.arange(0, 365, 5), ts['P']),
                      'diat': np.interp(np.arange(365), np.arange(0, 365, 5), ts['diat']),
                      'flag': np.interp(np.arange(365), np.arange(0, 365, 5), ts['flag']),
                      'microz': np.interp(np.arange(365), np.arange(0, 365, 5), ts['microz']),
                      'mesoz': np.interp(np.arange(365), np.arange(0, 365, 5), ts['mesoz']),
                      'ialg': np.interp(np.arange(365), np.arange(0, 365, 5), ts['ialg'])}
                
                return ts
            
            else:
                return ts
                
            
        if year1 < year2:
            
            ts2 = {'t_yr': list(np.tile(range(0,365), (year2-year1+1))),
                  't': list(range(0, 365*(year2-year1+1))),
                  'year': list(range(year1, year2+1)),
                  'T0': [],
                  'Td': [],
                  'P': [],
                  'diat': [],
                  'flag': [],
                  'microz': [],
                  'mesoz': [],
                  'ialg': []}
            
            ts['t'].extend(list(range(0, 73*(year2-year1+1))))
            ts['year'] = list(range(year1, year2+1))
            
            for yr in range(year1, year2+1):
                
                BB_profiles = BB_profiles_all[yr]
                
                # Time
                ts['t_yr'].extend(list(range(0, BB_profiles['mldprof'].shape[1])))
                
                # Surface Temperature (° Celsius)
                temp = BB_profiles['TEMsurf']
                depth = BB_profiles['depth']
                
                ts['T0'].extend(np.nanmean(np.where(depth <= 100, temp, np.nan), axis=0))
                
                # Deep Temperature (° Celsius)
                ts['Td'].extend(np.nanmean(np.where(depth >= 300, temp, np.nan), axis=0))
                
                # Prey (mg.chl.m-3)
                prey = BB_profiles['dia2prof'] + BB_profiles['fla2prof']
                
                ts['P'].extend(np.nanmax(prey, axis=0))
                
                # Diatoms (mg.chl.m-3)
                diat = BB_profiles['dia2prof']
                ts['diat'].extend(np.nanmax(diat, axis=0))
                
                # Flagellates (mg.chl.m-3)
                flag = BB_profiles['fla2prof']
                ts['flag'].extend(np.nanmax(flag, axis=0))
                
                # Microzooplankton (mmolN.m-3)
                microz = BB_profiles['micprof']
                ts['microz'].extend(np.nanmean(np.where(depth <= 100, microz, np.nan), axis=0))
                
                # Mesozooplankton (mmolN.m-3)
                mesoz = BB_profiles['mesprof']
                ts['mesoz'].extend(np.nanmean(np.where(depth <= 100, mesoz, np.nan), axis=0))
                
                # Ice Algae (mmolN.m-3)
                ialg = BB_profiles['ialgprof']
                ts['ialg'].extend(np.nanmax(ialg, axis=0))
                
                if extend:
                
                    initial_idx = np.arange(0, 365, 5)
                    new_idx = np.arange(365)
                    
                    ts2['T0'].extend(np.interp(new_idx, initial_idx, ts['T0'][-73:]))
                    ts2['Td'].extend(np.interp(new_idx, initial_idx, ts['Td'][-73:]))
                    ts2['P'].extend(np.interp(new_idx, initial_idx, ts['P'][-73:]))
                    ts2['diat'].extend(np.interp(new_idx, initial_idx, ts['diat'][-73:]))
                    ts2['flag'].extend(np.interp(new_idx, initial_idx, ts['flag'][-73:]))
                    ts2['microz'].extend(np.interp(new_idx, initial_idx, ts['microz'][-73:]))
                    ts2['mesoz'].extend(np.interp(new_idx, initial_idx, ts['mesoz'][-73:]))
                    ts2['ialg'].extend(np.interp(new_idx, initial_idx, ts['ialg'][-73:]))
                    
            if extend:
                return ts2
        
            else:
                return ts
                
                    # ts = {'t_yr': list(np.tile(range(0,365), (year2-year1+1))),
                    #       't': list(range(0, 365*(year2-year1+1))),
                    #       'year': list(range(year1, year2+1)),
                    #       'T0': np.interp(new_idx, initial_idx, ts['T0']),
                    #       'Td': np.interp(new_idx, initial_idx, ts['Td']),
                    #       'P': np.interp(new_idx, initial_idx, ts['P']),
                    #       'diat': np.interp(new_idx, initial_idx, ts['diat']),
                    #       'flag': np.interp(new_idx, initial_idx, ts['flag']),
                    #       'microz': np.interp(new_idx, initial_idx, ts['microz']),
                    #       'mesoz': np.interp(new_idx, initial_idx, ts['mesoz']),
                    #       'ialg': np.interp(new_idx, initial_idx, ts['ialg'])}
    
    if region == 'NOW':
        
        # Open the dict with all the profiles (all years)
        script_dir = os.path.dirname(os.path.abspath(__file__))  # Répertoire du script
        file_path = os.path.join(script_dir, 'NOW_profiles_2002-2021_dict.pkl')
        
        with open(file_path, 'rb') as file:
                NOW_profiles_all = pickle.load(file)

        # Select only the years we want
        if year1 == year2:
            NOW_profiles = NOW_profiles_all[year1]
            
            # Time
            ts['t'] = list(range(0,365,5))
            ts['t_yr'] = list(range(0,365,5))
            ts['year'] = list(range(year1, year2+1))
            
            # Surface Temperature (° Celsius)
            temp = NOW_profiles['TEMsurf']
            depth = NOW_profiles['depth']
            
            ts['T0'] = np.nanmean(np.where(depth <= 100, temp, np.nan), axis=0)
            
            # Deep Temperature (° Celsius)
            ts['Td'] = np.nanmean(np.where(depth >= 300, temp, np.nan), axis=0)
            
            # Prey (mg.chl.m-3)
            prey = NOW_profiles['dia2prof'] + NOW_profiles['fla2prof']
            
            ts['P'] = np.nanmax(prey, axis=0)
            
            # Diatoms (mg.chl.m-3)
            diat = NOW_profiles['dia2prof']
            ts['diat'] = np.nanmax(diat, axis=0)
            
            # Flagellates (mg.chl.m-3)
            flag = NOW_profiles['fla2prof']
            ts['flag'] = np.nanmax(flag, axis=0)
            
            # Microzooplankton (mmolN.m-3)
            microz = NOW_profiles['micprof']
            ts['microz'] = np.nanmean(np.where(depth <= 100, microz, np.nan), axis=0)
            
            # Mesozooplankton (mmolN.m-3)
            mesoz = NOW_profiles['mesprof']
            ts['mesoz'] = np.nanmean(np.where(depth <= 100, mesoz, np.nan), axis=0)
            
            # Ice Algae (mmolN.m-3)
            ialg = NOW_profiles['ialgprof']
            ts['ialg'] = np.nanmax(ialg, axis=0)
            
            if extend:
                
                ts = {'t_yr': list(range(0,365)),
                      't': list(range(0,365)),
                      'year': list(range(year1, year2+1)),
                      'T0': np.interp(np.arange(365), np.arange(0, 365, 5), ts['T0']),
                      'Td': np.interp(np.arange(365), np.arange(0, 365, 5), ts['Td']),
                      'P': np.interp(np.arange(365), np.arange(0, 365, 5), ts['P']),
                      'diat': np.interp(np.arange(365), np.arange(0, 365, 5), ts['diat']),
                      'flag': np.interp(np.arange(365), np.arange(0, 365, 5), ts['flag']),
                      'microz': np.interp(np.arange(365), np.arange(0, 365, 5), ts['microz']),
                      'mesoz': np.interp(np.arange(365), np.arange(0, 365, 5), ts['mesoz']),
                      'ialg': np.interp(np.arange(365), np.arange(0, 365, 5), ts['ialg'])}
                
                return ts
            
            else:
                return ts
            
        if year1 < year2:
            
            ts2 = {'t_yr': list(np.tile(range(0,365), (year2-year1+1))),
                  't': list(range(0, 365*(year2-year1+1))),
                  'year': list(range(year1, year2+1)),
                  'T0': [],
                  'Td': [],
                  'P': [],
                  'diat': [],
                  'flag': [],
                  'microz': [],
                  'mesoz': [],
                  'ialg': []}
            
            ts['t'].extend(list(range(0, 73*(year2-year1+1))))
            ts['year'] = list(range(year1, year2+1))
            
            for yr in range(year1, year2+1):
                
                NOW_profiles = NOW_profiles_all[yr]
                
                # Time
                ts['t_yr'].extend(list(range(0, NOW_profiles['mldprof'].shape[1])))
                
                # Surface Temperature (° Celsius)
                temp = NOW_profiles['TEMsurf']
                depth = NOW_profiles['depth']
                
                ts['T0'].extend(np.nanmean(np.where(depth <= 100, temp, np.nan), axis=0))
                
                # Deep Temperature (° Celsius)
                ts['Td'].extend(np.nanmean(np.where(depth >= 300, temp, np.nan), axis=0))
                
                # Prey (mg.chl.m-3)
                prey = NOW_profiles['dia2prof'] + NOW_profiles['fla2prof']
                
                ts['P'].extend(np.nanmax(prey, axis=0))
                
                # Diatoms (mg.chl.m-3)
                diat = NOW_profiles['dia2prof']
                ts['diat'].extend(np.nanmax(diat, axis=0))
                
                # Flagellates (mg.chl.m-3)
                flag = NOW_profiles['fla2prof']
                ts['flag'].extend(np.nanmax(flag, axis=0))
                
                # Microzooplankton (mmolN.m-3)
                microz = NOW_profiles['micprof']
                ts['microz'].extend(np.nanmean(np.where(depth <= 100, microz, np.nan), axis=0))
                
                # Mesozooplankton (mmolN.m-3)
                mesoz = NOW_profiles['mesprof']
                ts['mesoz'].extend(np.nanmean(np.where(depth <= 100, mesoz, np.nan), axis=0))
                
                # Ice Algae (mmolN.m-3)
                ialg = NOW_profiles['ialgprof']
                ts['ialg'].extend(np.nanmax(ialg, axis=0))
                
                if extend:
                
                    initial_idx = np.arange(0, 365, 5)
                    new_idx = np.arange(365)
                    
                    ts2['T0'].extend(np.interp(new_idx, initial_idx, ts['T0'][-73:]))
                    ts2['Td'].extend(np.interp(new_idx, initial_idx, ts['Td'][-73:]))
                    ts2['P'].extend(np.interp(new_idx, initial_idx, ts['P'][-73:]))
                    ts2['diat'].extend(np.interp(new_idx, initial_idx, ts['diat'][-73:]))
                    ts2['flag'].extend(np.interp(new_idx, initial_idx, ts['flag'][-73:]))
                    ts2['microz'].extend(np.interp(new_idx, initial_idx, ts['microz'][-73:]))
                    ts2['mesoz'].extend(np.interp(new_idx, initial_idx, ts['mesoz'][-73:]))
                    ts2['ialg'].extend(np.interp(new_idx, initial_idx, ts['ialg'][-73:]))
                    
            if extend:
                return ts2
        
            else:
                return ts



