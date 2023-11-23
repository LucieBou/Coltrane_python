# -*- coding: utf-8 -*-

'''
Coltrane - yearday

@author Lucie Bourreau & Neil Banas
@date 2023/09/24
'''

def yearday(t):
    '''
    Give the day of the year.
    
     Parameters
    ----------
    t: array
        Dates in type datetime64 (e.g. t = np.array(['2023-08-18', '2023-08-19'], dtype='datetime64')) 
        or directly days of the year, in that case the output will be the same as t.
    
    Returns
    -------
    yday: array
        Associated day of the year.
    '''

    yday = t%365+1 # No sure the +1 is needed...
    
    return yday