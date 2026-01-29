# -*- coding: utf-8 -*-

'''
Coltrane - Variable Activity level example

@author Lucie Bourreau
@date 2025/10/07
'''

# Load the packages

import numpy as np
import matplotlib.pyplot as plt
import sys
import pickle
import pandas as pd

# Load the functions

sys.path.append("/Users/luciebourreau/Library/CloudStorage/OneDrive-UniversitéLaval/PhD_ULaval/Github_Lucie/Coltrane_python")

from coltrane_forcing import coltrane_forcing
from coltrane_params import coltrane_params
from coltrane_population import coltrane_population
from coltrane_community import coltrane_community

# Forcing
forcing = coltrane_forcing("NOW", 5)

# Parameters
#dt=30
dt=20
param = coltrane_params(requireActiveSpawning = 0,
                    a_winter = 0.2,
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
                    maxReserveFrac = 0.8509,
                    Ddia = 0)

#### Run coltrane_integrate line by line
# t0, s = timing_combinations(forcing, param)

# i = 0
# p = add_strategy_to_params(param, s, i)

### Run a population
pop, popts = coltrane_population(forcing, param, 2)

### Same but without active diapause

param2 = coltrane_params(requireActiveSpawning = 0,
                    a_winter = 0,
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
for ind in range(popts['D'].shape[1]):
    for strat in range(popts['D'].shape[2]):
        is_nan = np.isnan(popts['D'][:, ind, strat])
        mask = is_nan & np.cumsum(~is_nan).astype(bool)
        if np.any(mask):
            print(f"→ ind={ind}, strat={strat}, first NaN after non-NaN at t={np.where(mask)[0][0]}")
            break
    else:
        continue
    break

is_nan = np.isnan(popts['D'][:, ind, strat])
mask = is_nan & np.cumsum(~is_nan).astype(bool) # nan after it has been non-nan

is_nan2 = np.isnan(popts2['D'][:, ind, strat])
mask2 = is_nan2 & np.cumsum(~is_nan2).astype(bool) # nan after it has been non-nan

tdied_active = np.where(mask)[0][0] if np.any(mask) else None
tdied_diapause = np.where(mask2)[0][0] if np.any(mask2) else None


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

# Lipid reserves dans fullness

fig, ax = plt.subplots(1, 1, 
                       figsize=(12, 8))
    
# Generation Length
ax.scatter(popts['R'],
           popts['R']/popts['W'],
           label='Active diapause',
           edgecolors='#4c98fe',
           facecolors='#4c98fe',
           marker="o",
           alpha=0.5,
           s=np.nan_to_num(normalize(popts['dF1']),0) * 200)

ax.scatter(popts2['R'],
           popts2['R']/popts2['W'],
           label='Real diapause',
           edgecolors='#7fcdbb',
           facecolors='#7fcdbb',
           marker="o",
           alpha=0.5,
           s=np.nan_to_num(normalize(popts2['dF1']),0) * 200)

ax.set_ylabel(r'Lipid reserves ($\mu$gC)')
ax.set_xlabel('Lipid fullness')

ax.spines['top'].set_visible(False) 
ax.spines['right'].set_visible(False) 

plt.legend()
plt.show()

# Test forcing Qik 2015

forcing_qik = coltrane_forcing("Qik_mod_2015", 1)

months_labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

fig, ax = plt.subplots(2, 1, figsize=(8, 4), sharex=True)

ax[0].plot(forcing_qik['t'], forcing_qik['T0'], '-', c='skyblue', label='Surface temperature')
ax[0].plot(forcing_qik['t'], forcing_qik['Td'], '-', c='royalblue', label='Deep temperature')
ax[0].set_ylabel('Temperature (°C)', color='k')
ax[0].tick_params(axis='y', labelcolor='k')
ax[0].set_xticks(forcing_qik['t'][::31])
ax[0].set_xticklabels(months_labels)
ax[0].spines['top'].set_visible(False) 
ax[0].spines['right'].set_visible(False) 

ax[1].plot(forcing_qik['t'], forcing_qik['P'], '-', color='darkgreen', label='Chlorophylle', zorder=3)
ax[1].set_ylabel('Prey (mg.chl.m$^{-3}$)', color='k')
ax[1].tick_params(axis='y', labelcolor='k')
ax[1].spines['top'].set_visible(False) 
ax[1].spines['right'].set_visible(False) 

fig.tight_layout()
fig.legend(loc="center right", frameon=False, bbox_to_anchor=(1.25, 0.5))
plt.show()

# Test forcing Qik 2015 v2

forcing_qik2 = coltrane_forcing("Qik_obs_2015", 1)

fig, ax = plt.subplots(2, 1, figsize=(8, 4), sharex=True)

ax[0].plot(forcing_qik2['t'], forcing_qik2['T0'], '-', c='skyblue', label='Surface temperature')
ax[0].plot(forcing_qik2['t'], forcing_qik2['Td'], '-', c='royalblue', label='Deep temperature')
ax[0].set_ylabel('Temperature (°C)', color='k')
ax[0].tick_params(axis='y', labelcolor='k')
ax[0].set_xticks(forcing_qik2['t'][::31])
ax[0].set_xticklabels(months_labels)
ax[0].spines['top'].set_visible(False) 
ax[0].spines['right'].set_visible(False) 
#ax[0].legend(loc="center right", frameon=False, bbox_to_anchor=(1.5, 0.5))

ax[1].plot(forcing_qik2['t'], forcing_qik2['P'], '-', color='darkgreen', label='Chlorophylle', zorder=3)
ax[1].set_ylabel('Prey (mg.chl.m$^{-3}$)', color='k')
ax[1].tick_params(axis='y', labelcolor='k')
ax[1].spines['top'].set_visible(False) 
ax[1].spines['right'].set_visible(False) 
#ax[1].legend(loc="center right", frameon=False, bbox_to_anchor=(1.5, 0.5))

fig.tight_layout()
plt.show()


#### Figure multiple a_winter test

param0 = param.copy()
param0['a_winter'] = 0

param01 = param.copy()
param01['a_winter'] = 0.1

param02 = param.copy()

param03 = param.copy()
param03['a_winter'] = 0.3

param04 = param.copy()
param04['a_winter'] = 0.4

param05 = param.copy()
param05['a_winter'] = 0.5

param06 = param.copy()
param06['a_winter'] = 0.6

param1 = param.copy()
param1['a_winter'] = 1

pop0 = coltrane_population(forcing, param0, 1)
pop01 = coltrane_population(forcing, param01, 1)
pop02 = coltrane_population(forcing, param02, 1)
pop03 = coltrane_population(forcing, param03, 1)
pop04 = coltrane_population(forcing, param04, 1)
pop05 = coltrane_population(forcing, param05, 1)
pop06 = coltrane_population(forcing, param06, 1)
pop1 = coltrane_population(forcing, param1, 1)



fig, ax = plt.subplots(2, 2, 
                       figsize=(12, 8),
                       gridspec_kw={'wspace': 0.1, 'hspace': 0.15},
                       sharex='col')

ax = ax.flatten()
    
# Generation Length
ax[0].scatter(pop0['Wa'],
            (pop0['tEcen']-pop0['t0'])/365,
            label='Diapause',
            edgecolors='#969696',
            facecolors='#969696',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop0['F2']),0) * 200)
    
ax[2].scatter(pop01['Wa'],
            (pop01['tEcen']-pop01['t0'])/365,
            label='Activity = 0.1',
            edgecolors='#66c2a5',
            facecolors='#66c2a5',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop01['F2']),0) * 200)

ax[2].scatter(pop02['Wa'],
            (pop02['tEcen']-pop02['t0'])/365,
            label='Activity = 0.2',
            edgecolors='#fc8d62',
            facecolors='#fc8d62',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop02['F2']),0) * 200)

ax[2].scatter(pop03['Wa'],
            (pop03['tEcen']-pop03['t0'])/365,
            label='Activity = 0.3',
            edgecolors='#8da0cb',
            facecolors='#8da0cb',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop03['F2']),0) * 200)

ax[2].scatter(pop04['Wa'],
            (pop04['tEcen']-pop04['t0'])/365,
            label='Activity = 0.4',
            edgecolors='#e78ac3',
            facecolors='#e78ac3',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop04['F2']),0) * 200)

ax[2].scatter(pop05['Wa'],
            (pop05['tEcen']-pop05['t0'])/365,
            label='Activity = 0.5',
            edgecolors='#a6d854',
            facecolors='#a6d854',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop05['F2']),0) * 200)

ax[2].scatter(pop06['Wa'],
            (pop06['tEcen']-pop06['t0'])/365,
            label='Activity = 0.6',
            edgecolors='#ffd92f',
            facecolors='#ffd92f',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop06['F2']),0) * 200)

ax[2].scatter(pop1['Wa'],
            (pop1['tEcen']-pop1['t0'])/365,
            label='Activity = 1',
            edgecolors='#e5c494',
            facecolors='#e5c494',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop1['F2']),0) * 200)

# Capital fraction of egg production
ax[1].scatter(pop0['Wa'],
            pop0['capfrac'],
            # label='Diapause',
            edgecolors='#969696',
            facecolors='#969696',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop0['F2']),0) * 200)
    
ax[3].scatter(pop01['Wa'],
            pop01['capfrac'],
            # label='Activity = 0.1',
            edgecolors='#66c2a5',
            facecolors='#66c2a5',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop01['F2']),0) * 200)

ax[3].scatter(pop02['Wa'],
            pop02['capfrac'],
            # label='Activity = 0.2',
            edgecolors='#fc8d62',
            facecolors='#fc8d62',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop02['F2']),0) * 200)

ax[3].scatter(pop03['Wa'],
            pop03['capfrac'],
            # label='Activity = 0.3',
            edgecolors='#8da0cb',
            facecolors='#8da0cb',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop03['F2']),0) * 200)

ax[3].scatter(pop04['Wa'],
            pop04['capfrac'],
            # label='Activity = 0.4',
            edgecolors='#e78ac3',
            facecolors='#e78ac3',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop04['F2']),0) * 200)

ax[3].scatter(pop05['Wa'],
            pop05['capfrac'],
            # label='Activity = 0.5',
            edgecolors='#a6d854',
            facecolors='#a6d854',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop05['F2']),0) * 200)

ax[3].scatter(pop06['Wa'],
            pop06['capfrac'],
            # label='Activity = 0.5',
            edgecolors='#ffd92f',
            facecolors='#ffd92f',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop06['F2']),0) * 200)

ax[3].scatter(pop1['Wa'],
            pop1['capfrac'],
            # label='Activity = 0.5',
            edgecolors='#e5c494',
            facecolors='#e5c494',
            marker="o",
            alpha=0.3,
            s=np.nan_to_num(normalize(pop1['F2']),0) * 200)

max_value = max(np.nanmax(pop0['Wa']), 
                np.nanmax(pop01['Wa']), 
                np.nanmax(pop02['Wa']), 
                np.nanmax(pop03['Wa']), 
                np.nanmax(pop04['Wa']), 
                np.nanmax(pop05['Wa']),
                np.nanmax(pop06['Wa']),
                np.nanmax(pop1['Wa']))

ax[0].set_ylabel('Generation length (years)')
ax[0].set_ylim([0,4])

ax[2].set_ylabel('Generation length (years)')
ax[2].set_ylim([0,4])

ax[1].set_ylabel(r'Capital fraction of $E$')
ax[1].set_ylim([0,1])

ax[3].set_ylabel(r'Capital fraction of $E$')
ax[3].set_ylim([0,1])

ax[2].set_xlabel(r'Adult Body Size ($\mu$gC)')
ax[3].set_xlabel(r'Adult Body Size ($\mu$gC)')
ax[2].set_xlim([-0.2,max_value+100])
ax[3].set_xlim([-0.2,max_value+100])

ax[0].spines['top'].set_visible(False) 
ax[0].spines['right'].set_visible(False) 

ax[1].spines['top'].set_visible(False) 
ax[1].spines['right'].set_visible(False)


ax[2].spines['top'].set_visible(False) 
ax[2].spines['right'].set_visible(False)


ax[3].spines['top'].set_visible(False) 
ax[3].spines['right'].set_visible(False) 


fig.legend(
    loc='center right',
    bbox_to_anchor=(1.04, 0.5),
    frameon=False,              
)
# plt.tight_layout(rect=[0, 0, 0.9, 1])

plt.show()


################################################### Test run coltrane_community


traits = {}
#traits['a_winter'] = [0, 0.05, 0.1, 0.5]
traits['a_winter'] = [0, 0.5]

coltrane_community('winter_act_ex_test', forcing, param, traits)

path = '/Users/luciebourreau/Library/CloudStorage/OneDrive-UniversitéLaval/PhD_ULaval/Github_Lucie/Coltrane_python/winter_act_ex_test'
path = '/Users/luciebourreau/winter_act_ex_test'

with open(path, 'rb') as file:
    loaded_data = pickle.load(file)

comm = loaded_data['comm']


## F2 en fonction de tEcen 

tEcen = comm['tEcen']
tEcen_year = tEcen % 365
tEcen_year[tEcen_year == 0] = 365

F2 = comm['F2']

a_winter_vals = traits['a_winter']
#colors = ['#80b1d3', '#8dd3c7', '#bebada', '#fb8072']

colors = ['#80b1d3', '#fb8072']


fig, ax = plt.subplots(figsize=(8, 5))

for i, (aw, col) in enumerate(zip(a_winter_vals, colors)):
    ax.scatter(
        tEcen_year[i].ravel(),
        F2[i].ravel(),
        color=col,
        alpha=0.4,
        s=8,
        label=f"a_winter = {aw}"
    )
ax.set_xlabel("tEcen")
ax.set_ylabel("Fitness (F2)")
ax.legend(frameon=False, loc='upper right')
plt.show()


fig, ax = plt.subplots(figsize=(8, 5))

for i, (aw, col) in enumerate(zip(a_winter_vals, colors)):
    # Aplatir pour avoir tous les points de ce a_winter
    t_flat = tEcen_year[i].ravel()
    F_flat = F2[i].ravel()
    
    # Regrouper par jour de l’année et calculer la moyenne
    df = pd.DataFrame({'day': t_flat, 'fitness': F_flat})
    mean_by_day = df.groupby('day', as_index=False)['fitness'].mean()

    # Appliquer un lissage (rolling mean sur ±5 jours)
    mean_by_day['fitness_smooth'] = (
        mean_by_day['fitness']
        .rolling(window=40, center=True, min_periods=1)
        .mean()
    )

    # Tracer la courbe lissée
    ax.plot(mean_by_day['day'], mean_by_day['fitness_smooth'],
            color=col, lw=2, label=f"a_winter = {aw}")

ax.set_xlabel("tEcen")
ax.set_ylabel("Fitness (F2)")
ax.legend(frameon=False, loc='upper right')
plt.show()


## a_winter en fonction de F2


fig, ax = plt.subplots(figsize=(8, 5))

jitter_strength = 0.01

for i, (aw, col) in enumerate(zip(a_winter_vals, colors)):

    awinter_flat = comm['a_winter'][i].ravel()
    F2_flat = F2[i].ravel()

    tdia_enter = comm['tdia_enter'][i].ravel()
    tdia_exit  = comm['tdia_exit'][i].ravel()

    # Diapause duration (days)
    tdia_duration = np.where(tdia_exit >= tdia_enter,
                             tdia_exit - tdia_enter,
                             365 - tdia_enter + tdia_exit)
    
    # Jitter on winter activity level
    awinter_jittered = awinter_flat + np.random.uniform(
        low=-jitter_strength, high=jitter_strength, size=awinter_flat.shape
    )

    sc = ax.scatter(F2_flat, awinter_jittered,
                    c=tdia_duration,
                    cmap='viridis',
                    s=8,
                    lw=0.5,
                    alpha=0.8,
                    edgecolor='none')


ax.set_xlabel("Fitness (F2)")
ax.set_ylabel("Winter Activity Level")
cbar = plt.colorbar(sc, ax=ax, label="Diapause duration (days)")
plt.show()


## Durée diapause en fonction de F2

fig, ax = plt.subplots(figsize=(8, 5))

jitter_strength = 3

for i, (aw, col) in enumerate(zip(a_winter_vals, colors)):

    awinter_flat = comm['a_winter'][i].ravel()
    F2_flat = F2[i].ravel()

    tdia_enter = comm['tdia_enter'][i].ravel()
    tdia_exit  = comm['tdia_exit'][i].ravel()

    # Diapause duration (days)
    tdia_duration = np.where(tdia_exit >= tdia_enter,
                             tdia_exit - tdia_enter,
                             365 - tdia_enter + tdia_exit)
    
    # Jitter on winter activity level
    diaduration_jittered = tdia_duration + np.random.uniform(
        low=-jitter_strength, high=jitter_strength, size=tdia_duration.shape
    )

    ax.scatter(diaduration_jittered, F2_flat,
                    c=col,
                    #cmap='viridis',
                    s=8,
                    lw=0.5,
                    alpha=0.8,
                    edgecolor='none',
                    label=f"a_winter = {aw}")

ax.set_xlabel("Diapause duration (days)")
ax.set_ylabel("Fitness (F2)")
ax.legend(frameon=False, loc='upper left', markerscale=2, fontsize = 8)
plt.show()


# Genlen et capfrac


rows = []
for i, a_w in enumerate(a_winter_vals):
 
        tEcen = comm['tEcen'][i, :, :]
        t0 = comm['t0'][i, :, :]
        capfrac = comm['capfrac'][i, :, :]
        F2 = comm['F2'][i, :, :]
        Wa = comm['Wa'][i, :, :] 
    
        rows.append(pd.DataFrame({
            'a_winter': a_w,
            'Wa': Wa.ravel(),
            'genlen': ((tEcen - t0) / 365).ravel(),
            'capfrac': capfrac.ravel(),
            'F2': F2.ravel()
        }))
    
df = pd.concat(rows, ignore_index=True)
df = df.dropna(subset=['Wa', 'genlen', 'capfrac', 'F2'])
df = df[df['F2'] > 0]

df['size'] = 2 + 180 * (df['F2'] / df['F2'].max())
df['alpha'] = 0.3 + 0.7 * (df['F2'] / df['F2'].max())


fig, ax = plt.subplots(2, 2, 
                       figsize=(12, 8),
                       gridspec_kw={'wspace': 0.1, 'hspace': 0.15},
                       sharex='col')

ax = ax.flatten()
    
# Generation Length
ax[0].scatter(df[df['a_winter'] == 0]['Wa'],
            df[df['a_winter'] == 0]['genlen'],
            label='Diapause',
            edgecolors='#969696',
            facecolors='#969696',
            marker="o",
            alpha=0.3,
            s=df[df['a_winter'] == 0]['size'])

ax[2].scatter(df[df['a_winter'] == 0.5]['Wa'],
            df[df['a_winter'] == 0.5]['genlen'],
            label='Activity = 0.5',
            edgecolors='#a6d854',
            facecolors='#a6d854',
            marker="o",
            alpha=0.3,
            s=df[df['a_winter'] == 0.5]['size'])

ax[1].scatter(df[df['a_winter'] == 0]['Wa'],
            df[df['a_winter'] == 0]['capfrac'],
            label='Diapause',
            edgecolors='#969696',
            facecolors='#969696',
            marker="o",
            alpha=0.3,
            s=df[df['a_winter'] == 0]['size'])

ax[3].scatter(df[df['a_winter'] == 0.5]['Wa'],
            df[df['a_winter'] == 0.5]['capfrac'],
            label='Activity = 0.5',
            edgecolors='#a6d854',
            facecolors='#a6d854',
            marker="o",
            alpha=0.3,
            s=df[df['a_winter'] == 0.5]['size'])


ax[0].set_ylabel('Generation length (years)')
ax[0].set_ylim([0,4])

ax[2].set_ylabel('Generation length (years)')
ax[2].set_ylim([0,4])

ax[1].set_ylabel(r'Capital fraction of $E$')
ax[1].set_ylim([0,1])

ax[3].set_ylabel(r'Capital fraction of $E$')
ax[3].set_ylim([0,1])

ax[2].set_xlabel(r'Adult Body Size ($\mu$gC)')
ax[3].set_xlabel(r'Adult Body Size ($\mu$gC)')
ax[2].set_xlim([-0.2,4000+100])
ax[3].set_xlim([-0.2,4000+100])

ax[0].spines['top'].set_visible(False) 
ax[0].spines['right'].set_visible(False) 

ax[1].spines['top'].set_visible(False) 
ax[1].spines['right'].set_visible(False)


ax[2].spines['top'].set_visible(False) 
ax[2].spines['right'].set_visible(False)


ax[3].spines['top'].set_visible(False) 
ax[3].spines['right'].set_visible(False) 


fig.legend(
    loc='center right',
    bbox_to_anchor=(1.04, 0.5),
    frameon=False,              
)
# plt.tight_layout(rect=[0, 0, 0.9, 1])

plt.show()
