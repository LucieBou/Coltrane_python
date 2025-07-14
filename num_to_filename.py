# -*- coding: utf-8 -*-

'''
Coltrane - num_to_filename

@author Lucie Bourreau & Neil Banas
@date 2023/08/24
'''

import os

def num_to_filename(basedir, n):
    '''
    Convention for converting run numbers to filenames in a nested structure that keeps directories 
    from filling up with tens of thousands of files.
    If n has more than 3 numbers the other numbers will create a new folder
    
    Example: 
    num2filename('my/dir/',345) gives my/dir/0/345.txt
    num2filename('my/dir/',12345) gives my/dir/12/345.txt
    
    Parameters
    ----------
    basedir: directory
        Directory to store the file.
    n: int
        Number to give to the filename.
    
    Returns
    -------
    fname: directory
        Directory with the name of the file (in .txt) depending on n.
    
    '''
    
    nstr = str(n)
    if len(nstr) < 4:
        nstr = '0000' + nstr
        nstr = nstr[-4:]
    
    fname = os.path.join(basedir, nstr[:-3], nstr[-3:] + '.txt')
    
    return fname
