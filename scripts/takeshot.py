#! /usr/bin/env python
# takeshot.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Take a 'shot' in FFS parlance, from one interface to another.  See
Allen, Valerani, ten Wolde J. Phys. Condens. matter 21, 463102.
"""

import sys
import os
import numpy as np
import funcselector
import readwrite
from ffsfunctions import *

# get command line arguments and complain if not right. We expect at
# least 1 argument: intfrom, which is the interface we are arriving
# from e.g. if intfrom=0, we are taking a shot from lambda0, aiming
# for lambda1. The second argument is optional, if present check if we
# already have minsuccess successful shots at this interface, and
# terminate immediately if so (the latter functionality is provided so
# that we can kill an entire 'batch' of shots, so codeffs.py for more
# details).

argc = len(sys.argv)
if (argc != 2 and argc != 3):
    sys.exit('Error: takeshot.py expected either one or two arguments')
# interface we are coming from
intfrom = int(sys.argv[1])

# read general simulation parameters from file
params = getpickparams()

# if 3 args (the third arg can be anything), check if we've already
# had a sufficient number of 'successful' shots, and terminate if so.
if argc == 3:
    nsuccess = getnumsuccess(intfrom+1)
    if nsuccess >= params['minsuccess']:
        sys.exit('Already have {0} successful shots at this '
                 'interface'.format(params['minsuccess']))

# the environment variable 'SGE_TASK_ID' stores my shot number.  See
# the (chunky) Sun Grid Engine manual for more info!
myjobnm = int(os.environ['SGE_TASK_ID'])

# get shot dictionary from previous interface.  The shot dictionary
# stores the number of successful shots and their numbers (and some
# other stuff), which allows us to pick an initial configuration for
# the current shot.
shotdict = getshotdict(intfrom)

# pick an initial configuration at random from previous interface.
# First pick random number between 0 and total weight.
r = shotdict['nsuccesseff']*np.random.rand()
wcounter = 0.0
for (num, w) in zip(shotdict['successnumbers'],
                    shotdict['successweights']):
    wcounter = wcounter + w
    if (r <= wcounter):
        initnum = num
        break
initfile = 'pos{0}_{1}.xyz'.format(intfrom, initnum)

# print some diagnostic information handy for debugging
print ('I found {0} successful shots at the previous interface - you should'
       'check this is correct'.format(shotdict['nsuccess']))
print 'These are runs {0}'\
      .format(','.join([str(i) for i in shotdict['successnumbers']]))
print 'I have chosen the initial config {0}'.format(initfile)
params['restartfile'] = initfile

# take the shot (see ffsfunctions.py)
success, weight, time, positions = takeshot(initfile, intfrom, params)

# print out whether success/fail and time
if success:
    # we reached the next interface
    sucstring = 'SUCCESS'
else:
    sucstring = 'FAIL'
    # set weight in case we survived pruning attempts and still failed
    weight = 0
    
print 'Shot number {0} finished in time {1} with status {2}'\
      .format(myjobnm, time, sucstring)

# if I was successful, I need to save my config
if success:
    # get the desired writexyz function.
    writexyzfunc = funcselector.FuncSelector(params).WriteXyzFunc()
    writexyzfunc('pos{0}_{1}.xyz'.format(intfrom + 1, myjobnm),
                 positions, params)

# finally (whether success or fail), write the to shotsi_j.out.  These
# files are read by 'finish.py'. Here i is the interface we are trying
# to reach, and j is the shot number.  We write 3 numbers on single
# line: initialconfignumber timetaken success
fname = 'shots{0}_{1}.out'.format(intfrom + 1, myjobnm)
fout = open(fname, 'w')
fout.write('from time success weight\n{0} {1} {2:d} {3:.6f}\n'\
           .format(initnum, time, success, weight))
fout.close()
