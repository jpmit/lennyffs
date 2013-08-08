# analyse.py
# James Mithen
# j.mithen@surrey.ac.uk

# functions for analysing FFS results

import os
import numpy as np

def ancestry(fname,names=True):
    """Return ancestry of FFS configuration"""
    posstr = fname.split('/')[-1]
    nums = posstr[3:-4]
    ns = nums.split('_')
    intnum = int(ns[0])
    shotnum = int(ns[1])

    # open each interface file in turn
    shotnums = np.zeros(intnum+1,dtype=int)
    shotnums[-1] = shotnum

    for i in range(intnum):
        fint = open('interface%d.out' %(intnum - i),'r')
        lines = fint.readlines()
        fint.close()
        last = shotnums[intnum - i]
        shotfrom = int(lines[last + 5].split()[1])
        shotnums[intnum-i-1] = shotfrom

    if not names:
        return shotnums[:-1]

    # stitch together filenames
    fnames = ['']*intnum
    for i in range(intnum):
        fnames[i] = 'pos%d_%d.xyz' %(i,shotnums[i])
    return fnames
