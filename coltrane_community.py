# -*- coding: utf-8 -*-

'''
Coltrane - coltrane_community

@author Lucie Bourreau & Neil Banas
@date 2023/11/03
'''

import numpy as np
import pickle

from coltrane_population import coltrane_population

def coltrane_community(outfile,forcing,p0,traits):
    '''
    Runs Coltrane for a single forcing time series but some range of one or more traits, 
    constructing a population for each element of the matched fields in traits. 
    
    Each population has a range of spawning dates t0 and timing strategy s = (tdia_enter, tdia_exit, dtegg).
    
    Parameters
    ----------
    outfile: string
        Name of the file to save.
    forcing: dict
        Set of forcing, composed of prey and temperature (surface & deep) cycle over several years.
    p0: dict
        Set of parameters (from coltrane_params.py)
    traits: dict
        Traits values to test.
    
    Returns
    -------
    None
    '''
    
    traitFields = list(traits.keys())
    Ntr = len(traits[traitFields[0]])

    print(f'Varying the traits {traitFields} in {Ntr} combinations.')

    # Initialize the output structure ----------------------------------------------------------------------        
    
    comm = {}

    # Run each trait combination ---------------------------------------------------------------------------       
    
    for i in range(Ntr):
        print('Trait combination', i + 1, '/', Ntr)
        
        # Parameter set for trait combination i
        p = p0.copy()
        for k in traitFields:
            p[k] = traits[k][i]

        # Run one population
        pop_i = coltrane_population(forcing, p, 1)
    
        # Add results to the main structure (comm)
        fields = pop_i.keys()
    
        for k in fields:
            if k not in pop_i:
                continue

            pop_i_var_k = pop_i[k].copy()
            pop_i_var_k = pop_i_var_k.reshape((1,) + pop_i_var_k.shape)

            if i == 0:  # Initialize the shape of the output structure comm
                shape = (Ntr,) + pop_i_var_k.shape[1:]
                comm[k] = np.empty(shape)
                comm[k][:] = np.nan

            comm[k][i, ...] = pop_i_var_k  # Save one slice of comm

    # Clean up output --------------------------------------------------------------------------------------               

    # Add traits to the main output structure

    for k in traitFields:
        reshaped_traits_k = np.array(traits[k]).reshape(-1, 1, 1)
        comm[k] = np.tile(reshaped_traits_k, (1,) + comm['F1'].shape[1:]) # of shape Ntr x NC x NS

    # Save the data

    with open(outfile, 'wb') as file:
        pickle.dump({'comm': comm, 'forcing': forcing, 'p0': p0, 'traits': traits}, file)
        
        # To open the data for later, here are the command lines:
            # with open('Example', 'rb') as file:
                # loaded_data = pickle.load(file)

            # comm2 = loaded_data['comm']
            # forcing2 = loaded_data['forcing']
            # p02 = loaded_data['p0']
            # traits2 = loaded_data['traits']
