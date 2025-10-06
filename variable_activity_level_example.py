# -*- coding: utf-8 -*-

'''
Coltrane - Variable Activity level example

@author Lucie Bourreau
@date 2025/07/18
'''

# Load the packages

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load the functions

from coltrane_integrate import coltrane_integrate
from coltrane_forcing import coltrane_forcing
from coltrane_params import coltrane_params
from timing_combinations import timing_combinations
from add_strategy_to_params import add_strategy_to_params
from coltrane_population import coltrane_population
from yearday import yearday

# Forcing
forcing = coltrane_forcing("NOW", 7)

# Parameters
dt=30
param = coltrane_params(requireActiveSpawning = 0,
                    AllowActiveDiapause = 1,
                    mortality_penalty = 0,
                    tdia_exit = range(30,130,dt),
                    tdia_enter = range(250,365,dt),
                    min_genlength_years = 1,
                    max_genlength_years = 3,
                    dt_spawn = dt,
                    preySatVersion = 'now_icealg',
                    I0 = 0.3339,
                    Ks = 1.1796,
                    KsIA = 0.2777,
                    u0 = 0.006,
                    rm = 0.0698,
                    maxReserveFrac = 0.8509)

#### Run coltrane_integrate line by line
# t0, s = timing_combinations(forcing, param)

# i = 0
# p = add_strategy_to_params(param, s, i)

### Run a population
pop, popts = coltrane_population(forcing, param, 2)

## Figures to explore
popts['yday'] = yearday(popts['t'])

in_diapause = (popts['yday'] >= pop['tdia_enter'][np.newaxis, :, :]) | \
              (popts['yday'] <= pop['tdia_exit'][np.newaxis, :, :])
              
active = (popts['a'] == 1)
alive = (popts['D'] != np.nan)

violations = active & in_diapause

time_idx, ind_idx, strat_idx = np.where(violations)
ind_ActiveDiapause = np.unique(ind_idx)

### Same but without active diapause

param2 = coltrane_params(requireActiveSpawning = 0,
                    AllowActiveDiapause = 0,
                    mortality_penalty = 0.1,
                    tdia_exit = range(30,130,dt),
                    tdia_enter = range(250,365,dt),
                    min_genlength_years = 1,
                    max_genlength_years = 3,
                    dt_spawn = dt,
                    preySatVersion = 'now_icealg',
                    I0 = 0.3339,
                    Ks = 1.1796,
                    KsIA = 0.2777,
                    u0 = 0.006,
                    rm = 0.0698,
                    maxReserveFrac = 0.8509)

pop2, popts2 = coltrane_population(forcing, param2, 2)

### Figures

## State variables time series for ind 0
ind = 3
strat = np.unique(strat_idx[ind_idx == ind])[0]

is_nan = np.isnan(popts['D'][:, ind, strat])
mask = is_nan & np.cumsum(~is_nan).astype(bool) # nan after it has been non-nan

is_nan2 = np.isnan(popts2['D'][:, ind, strat])
mask2 = is_nan2 & np.cumsum(~is_nan2).astype(bool) # nan after it has been non-nan

tdied_active = np.where(mask)[0][0]
tdied_diapause = np.where(mask2)[0][0]


fig, (ax_top, ax_bottom) = plt.subplots(
    2, 1, figsize=(8, 8), sharex=True, gridspec_kw={'height_ratios': [1, 2]}
)

# --- Partie du haut (nouvel axe Y indépendant)
ax_top.plot(popts2['t'][:tdied_diapause+20, ind, strat], popts2['sat'][:tdied_diapause+20, ind, strat], color="lightgray")
ax_top.plot(popts2['t'][:tdied_active+20, ind, strat], popts2['sat'][:tdied_active+20, ind, strat], color="lightgray")
ax_top.set_ylabel("Prey saturation")

# --- Partie du bas (ta figure actuelle avec 3 axes y)
ax1 = ax_bottom
ax1.plot(popts2['t'][:tdied_diapause+20, ind, strat], popts2['D'][:tdied_diapause+20, ind, strat], color="red", alpha=0.5)
ax1.plot(popts['t'][:tdied_active+20, ind, strat], popts['D'][:tdied_active+20, ind, strat], color="#ff7f7f", alpha=0.5)
ax1.set_ylabel("Development", color="red")
ax1.tick_params(axis="y", labelcolor="red")

ax2 = ax1.twinx()
ax2.plot(popts2['t'][:tdied_diapause+20, ind, strat], popts2['R'][:tdied_diapause+20, ind, strat], color="green", alpha=0.5)
ax2.plot(popts['t'][:tdied_active+20, ind, strat], popts['R'][:tdied_active+20, ind, strat], color="#7fbf7f", alpha=0.5)
ax2.set_ylabel("Reserves", color="green")
ax2.tick_params(axis="y", labelcolor="green")

ax3 = ax1.twinx()
ax3.spines['right'].set_position(('outward', 60))
ax3.plot(popts2['t'][:tdied_diapause+20, ind, strat], popts2['lnN'][:tdied_diapause+20, ind, strat], color="blue", alpha=0.5)
ax3.plot(popts['t'][:tdied_active+20, ind, strat], popts['lnN'][:tdied_active+20, ind, strat], color="#7f7fff", alpha=0.5)
ax3.set_ylabel("Survivorship", color="blue")
ax3.tick_params(axis="y", labelcolor="blue")

# Lignes verticales de référence
ax1.axvline(pop2['tdia_exit'][ind, strat], color="black", linestyle=":", alpha=0.3)
ax1.axvline(pop2['tdia_enter'][ind, strat], color="black", linestyle=":", alpha=0.3)
ax1.axvline((pop2['t0'][ind, strat] + pop2['dtegg'][ind, strat]), color="black", linestyle="--", alpha=0.3)

ax1.set_xlabel("Time")

plt.tight_layout()
plt.show()

## Life strategies

def normalize(values):
    values = np.array(values, dtype=np.float64)
    normed = np.zeros_like(values)
    mask = ~np.isnan(values)

    if np.any(mask):  # éviter division par 0
        vmin = np.min(values[mask])
        vmax = np.max(values[mask])
        if vmax - vmin != 0:
            normed[mask] = (values[mask] - vmin) / (vmax - vmin)
    return normed

fig, ax = plt.subplots(1, 2, 
                       figsize=(12, 8),
                       gridspec_kw={'wspace': 0.1, 'hspace': 0.15})
    
# Generation Length
ax[0].scatter(pop['Wa'],
            (pop['tEcen']-pop['t0'])/365,
            label='Active diapause',
            edgecolors='#4c98fe',
            facecolors='#4c98fe',
            marker="o",
            alpha=0.5,
            s=np.nan_to_num(normalize(pop['F2']),0) * 200)
    
ax[0].scatter(pop2['Wa'],
            (pop2['tEcen']-pop2['t0'])/365,
            label='Real diapause',
            edgecolors='#7fcdbb',
            facecolors='#7fcdbb',
            marker="o",
            alpha=0.5,
            s=np.nan_to_num(normalize(pop2['F2']),0) * 200)

# Capital fraction of egg production
ax[1].scatter(pop['Wa'],
            pop['capfrac'],
            label='Active diapause',
            edgecolors='#4c98fe',
            facecolors='#4c98fe',
            marker="o",
            alpha=0.5,
            s=np.nan_to_num(normalize(pop['F2']),0) * 200)
    
ax[1].scatter(pop2['Wa'],
            pop2['capfrac'],
            label='Real diapause',
            edgecolors='#7fcdbb',
            facecolors='#7fcdbb',
            marker="o",
            alpha=0.5,
            s=np.nan_to_num(normalize(pop2['F2']),0) * 200)

max_value = max(np.nanmax(pop['Wa']), np.nanmax(pop2['Wa']))

ax[0].set_ylabel('Generation length (years)')
ax[0].set_ylim([0,4])

ax[1].set_ylabel(r'Capital fraction of $E$')
ax[1].set_ylim([0,1])

ax[0].set_xlabel(r'Adult Body Size ($\mu$gC)')
ax[1].set_xlabel(r'Adult Body Size ($\mu$gC)')
ax[0].set_xlim([-0.2,max_value+100])
ax[1].set_xlim([-0.2,max_value+100])

ax[0].spines['top'].set_visible(False) 
ax[0].spines['right'].set_visible(False) 

ax[1].spines['top'].set_visible(False) 
ax[1].spines['right'].set_visible(False) 

plt.legend()
plt.show()

## Delay between spawning and diapause entry that lead to active diapause
t0 = []
delta = []
for i in ind_ActiveDiapause:
    for s in np.unique(strat_idx[ind_idx == i]):
        yday_t0 = pop['t0'][i,s] % 365
        d = pop['tdia_enter'][i,s] - yday_t0
        t0.append(yday_t0)
        delta.append(d)
       
inactive = (popts['a'] == 0)
alive = (popts['D'] != np.nan)        
diapausing = alive & in_diapause & inactive
mask_diap = diapausing.any(axis=0)

plt.scatter(pop['t0'] % 365, pop['tdia_enter'] - pop['t0'] % 365, s = 10,color = 'lightgray')
plt.scatter((pop['t0'] % 365)[mask_diap], pop['tdia_enter'][mask_diap] - (pop['t0'] % 365)[mask_diap], s = 20,color = '#7fcdbb')
plt.scatter(t0, delta, c="#4c98fe", s=50)
plt.xlabel("Spawning date (days)", fontsize = 14)
plt.ylabel(r"$\Delta t = t_{\mathrm{dia\_enter}} - t_{0}$ (days)", fontsize = 14)
plt.show()
        