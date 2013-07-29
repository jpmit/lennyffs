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

import potselector
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

# from parameters file, create PotSelector object.  This will handle
# correct selection of the underlying fortran functions correctly (the
# functions called depend on the potential, i.e. the value of
# params['potential']).
PotManager = potselector.PotSelector(params)

# compute initial energy
epot = PotManager.TotalEnergy(positions,params)

# perform MC simulation
starttime = time.time()
params['cycle'] = params['ncycle']

if params['mctype'] == 'npt':
    positions, epot = PotManager.CycleNPT(positions,params,epot)
elif params['mctype'] == 'nvt':
    positions, epot = PotManager.CycleNVT(positions,params,epot)

endtime = time.time()

# write final positions to file
writeoutput.writexyz('finalpositions.xyz',positions,params)

# write runtime to stderr
sys.stderr.write("runtime in s: %.3f\n" %(endtime - starttime))
