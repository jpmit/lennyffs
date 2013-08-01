#! /usr/bin/env python
# code.py
# James Mithen
# j.mithen@surrey.ac.uk
#
# MC code for investigating epitaxial nucleation on surfaces.  This
# code just performs a certain number of MC steps governed by the
# parameters in the 'in' file.  Note there is no 'Order Parameter'
# calculation done here.  Thus, this script is mainly useful for
# equilibration/quenching etc.  For FFS simulations, see codeffs.py

import funcselector
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

# from parameters file, create FuncSelector object.  This will handle
# correct selection of the underlying fortran functions correctly (the
# functions called depend on the potential, i.e. the value of
# params['potential'], and also on the type of MC cycle wanted, i.e.
# params['mctype']
funcman = funcselector.FuncSelector(params)
totalenergy = funcman.TotalEnergyFunc()
runcycle = funcman.MCCycleFunc()

# compute initial energy
epot = totalenergy(positions,params)

# perform MC simulation
starttime = time.time()
params['cycle'] = params['ncycle']

# run the MC cycles
positions, epot = runcycle(positions, params, epot)

endtime = time.time()

# write final positions to file
writeoutput.writexyz('finalpositions.xyz',positions,params)

# write runtime to stderr
sys.stderr.write("runtime in s: %.3f\n" %(endtime - starttime))
