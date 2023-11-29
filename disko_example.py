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
from coltrane_forcing import coltrane_forcing

# Forcing --------------------------------------------------------------------------------------------------       

forcing = coltrane_forcing("DiskoBay", 7)

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

