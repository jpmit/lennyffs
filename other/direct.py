#!/usr/bin/env python
# direct.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Script for direct estimation of nulceation rate i.e. not using FFS.
At the moment this only works for Monte-Carlo integrators, but
extending it to Molecular Dynamics integrators would be
straightforward.
"""

import sys

import initsim
import writeoutput
import energy
import funcselector
from copy import deepcopy
from ffsfunctions import savelambda0config

# if true, we write the order parameter each time it is measured to
# stderr, otherwise we write order parameter to the file 'opval.out'.
# Writing to stderr is useful if we are running multiple instances of
# this script.  In both cases, we append the time at which we reach
# phase B (i.e. the mean first passage time) to 'times.out'
WRITE_TO_STDERR = True

if WRITE_TO_STDERR:
    opfile = sys.stderr
else:
    opfile = open('opval.out', 'a')

if len(sys.argv) != 2:
    sys.exit("Syntax direct.py OP")

# order parameter in product state
try:
    prodop = int(sys.argv[1])
except ValueError:
    sys.exit("OP must be an integer")

# read input parameters and write to file
params = initsim.getparams()
writeoutput.writepickparams(params)
writeoutput.writeparams(params)

# From params dictionary create FuncSelector object.  This will handle
# correct selection of the underlying fortran/C++ functions correctly
# (the functions called depend on the potential, i.e. the value of
# params['potential'], and also on the type of MC cycle wanted, i.e.
# params['mctype'], and on the order parameter desired,
# params['orderparam'].
funcman = funcselector.FuncSelector(params)
# TODO: totalenergy should be funcman.TotalEnergyFunc() but this does
# not seem to be working (?).
totalenergy = energy.len_totalenergy
runcycle = funcman.MCCycleFunc()
orderp = funcman.OrderParamFunc()
writexyz = funcman.WriteXyzFunc()

# initialize positions
positions = initsim.initpositions(params)

# we store the initial positions so that when we reach the product
# state we can start all over again.
initpositions = deepcopy(positions)

# write initial positions to file if new simulation
if params['simulation'] == 'new':
    writeoutput.writexyz('initpositions.xyz', positions, params)

# compute initial energy
epot = totalenergy(positions, params)
epotinit = epot

# write initial order parameter
opval = orderp(positions, params)
opfile.write('{0} {1}\n'.format(0, opval))
opfile.flush()

lamsamp = int(params['lambdasamp'])
params['cycle'] = lamsamp
cyclesdone = 0
qhits = 0

# we go to the product state 'totalqhits' times
while (qhits < int(params['totalqhits'])):
    # do some cycles
    positions, epot = runcycle(positions, params, epot)
    cyclesdone += lamsamp
    # compute order parameter
    opval = orderp(positions, params)
    opfile.write('{0} {1}\n'.format(cyclesdone, opval))
    opfile.flush()
    # if we reached the product state, record config and go back to
    # start
    if (opval >= prodop):
        qhits = qhits + 1
        # hacky: pass qhits + 1 to ensure we are always appending to
        # the 'times.out' file.  This means we can have multiple
        # instances of this script running without overwriting
        # 'times.out'.
        savelambda0config(qhits + 1, cyclesdone, opval, positions,
                          params, writexyz)
        # reset cyclesdone and positions
        cyclesdone = 0
        positions = initpositions
        epot = epotinit

opfile.close()
