#! /usr/bin/env python
# takeshot.py
# James Mithen
# j.mithen@surrey.ac.uk
#
# Take a 'shot' in FFS parlance, from one interface to another.  See
# Allen, Valerani, ten Wolde J. Phys. Condens. matter 21, 463102.

import sys
import os
import readwrite
from writeoutput import writexyz
import numpy as np
from ffsfunctions import *

# get arguments and complain if not right. Expect at least 1 argument:
# intfrom, which is the interface we are arriving from e.g. if
# intfrom=0, we are taking a shot from lambda0, aiming for
# lambda1. The second argument is optional, if present check if we
# already have minsuccess successful shots at this interface, and
# terminate if so

argc = len(sys.argv)
if (argc != 2 and argc != 3):
    sys.exit('Error: takeshot.py expected either one or two arguments')
intfrom = int(sys.argv[1])

# read general simulation parameters from file
params = getpickparams()

if argc == 3:
    nsuccess = getnumsuccess(intfrom+1)
    if nsuccess >= params['minsuccess']:
        sys.exit('Already have %d successful shots at this '
                 'interface' %params['minsuccess'])

# the environment variable 'SGE_TASK_ID' stores my shot number
myjobnm = int(os.environ['SGE_TASK_ID'])

# get shot dictionary from previous interface
shotdict = getshotdict(intfrom)

# pick an initial configuration at random from previous interface
# first pick random number between 0 and total weight
r = shotdict['nsuccesseff']*np.random.rand()
wcounter = 0.0
for (num, w) in zip(shotdict['successnumbers'],shotdict['successweights']):
    wcounter = wcounter + w
    if (r <= wcounter):
        initnum = num
        break
initfile = 'pos%d_%d.xyz' %(intfrom,initnum)

# print some diagnostic information handy for debugging
print ('I found %d successful shots at the previous interface - you should '
       'check this is correct' %shotdict['nsuccess'])
print 'These are runs %s' %','.join([str(i) for i in
                                     shotdict['successnumbers']])
print 'I have chosen the initial config %s' %initfile
params['restartfile'] = initfile

# take the shot
success,weight,time,positions = takeshot(initfile,intfrom,params)

# print out whether success/fail and time
if success:
    # we reached the next interface
    sucstring = 'SUCCESS'
else:
    sucstring = 'FAIL'
    # set weight in case we survived pruning attempts and still failed
    weight = 0 
print 'Shot number %d finished in time %d with status %s' %(myjobnm,time,sucstring)

# if I was successful, I need to save my config
if success:
    writexyz('pos%d_%d.xyz' %(intfrom+1,myjobnm),positions,params)

# finally (whether success or fail), write the to shotsi_j.out
# these files are read by finish.py
# here i=interface we are trying to reach, and j is the shot number
# write 3 numbers on single line:
# initialconfignumber timetaken success

fname = 'shots%d_%d.out' %(intfrom+1,myjobnm)
fout = open(fname,'w')
fout.write('from time success weight\n%d %d %d %.6f\n' %(initnum,time,success,weight))
fout.close()
