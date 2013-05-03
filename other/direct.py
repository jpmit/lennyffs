#!/usr/bin/env python
# direct.py
# 19th September 2012
# James Mithen
#
# Script for direct estimation of nulceation rate
# i.e. not using FFS

import initsim
import writeoutput
import energy
import mccycle
import bops
import numpy as N
import time
from copy import deepcopy
from bops import bopxbulk
from ffsfunctions import savelambda0config

ncrit = 150 # critical nucleus value

# read input parameters and write to file
params = initsim.getparams()
writeoutput.writepickparams(params)
writeoutput.writeparams(params)

# initialize positions
positions = initsim.initpositions(params)
initpositions = deepcopy(positions)

# write initial positions to file if new simulation
if params['simulation'] == 'new':
    writeoutput.writexyz('initpositions.xyz',positions,params)

# compute initial energy
epot = energy.totalenergy(positions,params)
epotinit = epot

# open bopx file
fout = open('bopx.out','w')

# perform MC simulation
nxtal, bopx = bopxbulk(positions,params)
fout.write("BOPX: %d\n" %bopx)
fout.flush()
lamsamp = int(params['lambdasamp'])
params['cycle'] = lamsamp
ttot = 0
qhits = 0
while (qhits < int(params['totalqhits'])):
    # do some cycles
    positions, epot = mccycle.cycle(positions,params,epot)
    ttot = ttot + lamsamp
    # compute bopxbulk
    nxtal, bopx = bopxbulk(positions,params)
    fout.write("BOPX: %d\n" %bopx)
    fout.flush()
    # if we made a critical nucleus, record config and go back to start
    if (bopx > ncrit):
        qhits = qhits + 1
        savelambda0config(qhits,ttot,positions,params)
        # reset ttot and positions
        ttot = 0
        positions = initpositions
        epot = epotinit

fout.close()
