#!/usr/bin/env python
# multisave.py
# 6th July 2013
# James Mithen
#
# Script for running simulations and saving NSAVE times
# Monte-Carlo cycles

import initsim
import writeoutput
import energy
import mccycle
import time
import sys

# number of MC cycles to run between saving
# the total number of MC cycles we run will be
# NSAVE*params['ncycle']
NSAVE = 1000

# read input parameters and write to file
params = initsim.getparams()
writeoutput.writepickparams(params)
writeoutput.writeparams(params)

# initialize positions
positions = initsim.initpositions(params)

# write initial positions to file if new simulation
if params['simulation'] == 'new':
    writeoutput.writexyz('initpositions.xyz',positions,params)

# compute initial energy
epot = energy.totalenergy(positions,params)

# perform MC simulation
params['starttime'] = time.time()
params['cycle'] = params['ncycle']

if params['mctype'] == 'npt':
    cyclefunc = mccycle.cyclenpt
elif params['mctype'] == 'nvt':
    cyclefunc = mccycle.cycle

for rnum in range(1,NSAVE+1):
    positions, epot = cyclefunc(positions,params,epot)
    # write final positions to file
    writeoutput.writexyz('savepositions%d.xyz' %rnum,positions,params)
    # write params.pkl for every file
    # we do this since for the NPT simulations the box dimensions could
    # have changed (they almost certainly should have).
    writeoutput.writepickparams(params, 'params%d.pkl' %rnum)
    writeoutput.writeparams(params, 'params%d.out' %rnum)

params['endtime'] = time.time()

# write runtime to stderr
sys.stderr.write("runtime in s: %.3f" %(params['endtime']-
                                        params['starttime']))
