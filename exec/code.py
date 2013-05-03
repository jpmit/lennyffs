#! /usr/bin/env python
# code.py
# 26th June 2012
# James Mithen
#
# MC code for investigating epitaxial nucleation on surfaces
# this code just performs a certain number of MC steps governed
# by the parameters in the 'in' file.
# Note there is no 'Order Parameter' calculation done here
# Thus, this script is mainly useful for equilibration/quenching etc.
# For FFS simulations, see codeffs.py

import initsim
import writeoutput
import energy
import mccycle
import time
import sys

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
    positions, epot = mccycle.cyclenpt(positions,params,epot)
elif params['mctype'] == 'nvt':
    positions, epot = mccycle.cycle(positions,params,epot)

params['endtime'] = time.time()

# write final positions to file
writeoutput.writexyz('finalpositions.xyz',positions,params)

# write runtime to stderr
sys.stderr.write("runtime in s: %.3f" %(params['endtime']-params['starttime']))
