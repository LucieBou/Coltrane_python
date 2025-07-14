# Coltrane: from matlab to python
**Lucie Bourreau - November, 2023**

## Introduction

Here you can find all the functions needed to run the Coltrane 2.0 model in Python. Coltrane was originally coded in Matlab and whose original code can be found on Neil Banas' Github (https://github.com/neilbanas/coltrane). Coltrane was first published in 2016 (https://doi.org/10.3389/fmars.2016.00225) and has since received some improvements, which are explained in the headings of the relevant functions.   

Coltrane: Copepod Life-history Traits and Adaptation to Novel Environments.    

## How it works

Coltrane is an individual trait-based model that simulates life-history strategies of copepods in response to an environment characterized by water temperature and prey availability. It generates copepod-like agents (or cohorts) with a specific spawning date. Each cohort is assigned a suite of life strategies (diapause entry date, diapause exit date, date of egg production) and a set of state variables and timing metrics is computed to establish which life strategies are viable and which are optimal (regarding cohort's fitness), depending on the environment. The cohort’s life history can then be assessed by looking back to the trajectories of each traits/state variable.    

Each function is explained in the docstring but here is a quick overview:

*coltrane_integrate.py* - Heart of the model where the state variables and timing metrics are computed for each cohort at a specific life-strategy during all the running time we decide (for example, 7 years).        

*coltrane_population.py* - Run coltrane_integrate for each life-strategy to create a population and save the time series of state variables and timing metrics into two different outputs (time series are not always useful and are very time consuming).      

*coltrane_community.py* - Run coltrane_population for a range of traits values, for example by varying the value of the parameter "u0: food-saturated development rate at T = 0" and give the same outputs as coltrane_population but for every trait value.     

*coltrane_params.py* - Contains the default parameters of the model and can be use to make changes in the paramosome.       

*disko_example.py* - Example on how to use Coltrane with specific forcing and a paramosome. This example is explain in the paper of 2016 (https://doi.org/10.3389/fmars.2016.00225).   

**Please, do not hesitate to contact me at lucie.bourreau.1@ulaval.ca or Neil Banas at neil.banas@strath.ac.uk for any questions, comments or suggestions.**    
**You can also visit this website for more information about Coltrane: https://neilbanas.com/projects/coltrane/**

