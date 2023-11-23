# -*- coding: utf-8 -*-

'''
Coltrane - Disko Bay example (see Banas et al., 2016)

@author Lucie Bourreau & Neil Banas
@date 2023/11/20
'''

import numpy as np
import pickle
import matplotlib.pyplot as plt


from coltrane_community import coltrane_community
from coltrane_params import coltrane_params
from coltrane_population import coltrane_population

# Forcing --------------------------------------------------------------------------------------------------       

forcing = {}
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

## Repeat this cycle for 7 years to resolve 3 year life cycles

Nyears = 7
forcing['t'] = np.arange(0, Nyears * 365)
forcing['T0'] = np.tile(forcing['T0'], Nyears)
forcing['Td'] = np.tile(forcing['Td'], Nyears)
forcing['P'] = np.tile(forcing['P'], Nyears)

# Parameters and traits ------------------------------------------------------------------------------------        

dt = 20

p00 = coltrane_params(requireActiveSpawning = 0,
                    tdia_exit = range(80,120,dt),
                    tdia_enter = range(260,300,dt),
                    min_genlength_years = 1,
                    max_genlength_years = 3,
                    dt_spawn = 15,
                    preySatVersion = 'default',
                    I0 = 0.36,
                    Ks = 1,
                    maxReserveFrac = 0.8)

traits = {}
traits['u0'] = [0.005, 0.007, 0.009]

# Main Coltrane experiment --------------------------------------------------------------------------------- 

coltrane_community('disko_ex', forcing, p00, traits)

# Load outputs ---------------------------------------------------------------------------------------------

with open('disko_ex', 'rb') as file:
    loaded_data = pickle.load(file)

comm = loaded_data['comm']
p01 = loaded_data['p0']

comm['gl'] = (comm['tEcen'] - comm['t0']) / 365
comm['Fyr'] = comm['F2'] ** (1 / 2 / comm['gl'])

f = np.where((comm['Fyr'] > 0) & (comm['t0'] <= 365))

# Figures --------------------------------------------------------------------------------------------------

## Forcing

plt.figure(figsize=(12, 8))

plt.suptitle('Forcing', fontsize=16)

plt.subplot(311)
plt.scatter(forcing['t'] / 365, forcing['T0'], marker='o', s=10)
plt.xlabel('Time (year)')
plt.ylabel('Surface temperature (°C)')

plt.subplot(312)
plt.scatter(forcing['t'] / 365, forcing['Td'], marker='o', s=10)
plt.xlabel('Time (year)')
plt.ylabel('Deep temperature (°C)')

plt.subplot(313)
plt.scatter(forcing['t'] / 365, forcing['P'], marker='o', s=10)
plt.xlabel('Time (year)')
plt.ylabel('Prey (mg chl m-3)')

plt.tight_layout()
plt.show()

## Exploration

plt.figure(figsize=(12, 8))

plt.suptitle('Exploration', fontsize=16)

plt.subplot(221)
plt.scatter(comm['Wa'][f], comm['gl'][f], s=30, c=comm['Fyr'][f], marker='o', cmap='viridis')
plt.xscale('log')
plt.colorbar()
plt.xlabel('Wa')
plt.ylabel('gl')

plt.subplot(222)
plt.scatter(comm['Wa'][f], comm['Ra'][f] / comm['Wa'][f], s=30, c=comm['Fyr'][f], marker='o', cmap='viridis')
plt.xscale('log')
plt.colorbar()
plt.xlabel('Wa')
plt.ylabel('Ra/Wa')

plt.subplot(223)
plt.scatter(comm['Wa'][f], comm['D_winter'][f], s=30, c=comm['Fyr'][f], marker='o', cmap='viridis')
plt.xscale('log')
plt.colorbar()
plt.xlabel('Wa')
plt.ylabel('Dwinter')

plt.subplot(224)
plt.scatter(comm['Wa'][f], comm['capfrac'][f], s=30, c=comm['Fyr'][f], marker='o', cmap='viridis')
plt.xscale('log')
plt.colorbar()
plt.xlabel('Wa')
plt.ylabel('capfrac')

plt.tight_layout()
plt.show()

