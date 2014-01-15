#! /usr/bin/env python
# lambda0.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Take the system from phase A to Lambda0 in FFS simulation.  see Allen,
 Valerani, ten Wolde J. Phys. Condens. matter 21, 463102.
"""

import sys
import pickle
import numpy as np
import initsim
import energy
import mccycle
from ffsfunctions import savelambda0config, getpickparams
import potselector
from lenexceptions import *

params = getpickparams()

lamA = params['lambdaA']
lam0 = params['lambdas'][0]
totalqhits = params['totalqhits'] # num times to go through lambda0
                                  # from phase A
lamsamp = params['lambdasamp']    # num MC cycles per OP evaluation
params['cycle'] = lamsamp

# initialize positions
positions = initsim.initpositions(params)

# get the correct energy function and MC cycle function using
# PotSelector interface
potman = potselector.PotSelector(params)
totalenergyfunc = potman.TotalEnergyFunc()
cyclefunc = potman.MCCycleFunc()
orderpfunc = potman.OrderParamFunc()

# check that lambda < lamA (we are in phase A)
op = orderpfunc(positions,params)
if (op >= lamA):
    raise FFSError, ('OP is {0}, system must start in phase A, '
                     'OP < {1}'.format(op, lamA))

# write initial OP to file
fin = open('ffsacc.out', 'w')
fin.write('ncycles OP\n')
fin.write('{0} {1}\n'.format(0, op))
fin.flush()

# start simulation
qhits = 0 # num times passed through lambda0
thit = 0
ttot = 0

# initial potential energy
epot = totalenergyfunc(positions, params)

# evolve the system in time until it has 'hit' the interface (lambda0)
# the desired number of times.
while (qhits < totalqhits):
    # go forward in time
    # TODO: make this work for MD (will require velocities and forces)
    positions, epot = cyclefunc(positions, params, epot)
    thit = thit + lamsamp
    ttot = ttot + lamsamp
    # evaluate OP
    op = orderpfunc(positions, params)
    # write OP and total time to file
    fin.write('{0} {1} {2}\n'.format(ttot, nxtal, bopx))
    fin.flush()
    if (bopx >= lam0):
        # we hit the interface
        qhits = qhits + 1
        savelambda0config(qhits, thit, positions, params)
        thit = 0
        # now let the system relax back to lambda A
        while (bopx >= lamA):
            # TODO: make this work for MD (will require velocities and forces)
            positions, epot = cyclefunc(positions, params, epot)
            ttot = ttot + lamsamp                
            # evaluate OP
            op = orderpfunc(positions, params)
            # write OP and total time to file
            fin.write('{0} {1} {2}\n'.format(ttot, nxtal, bopx))
            # if we reached lam0*10 then return to phase A by
            # reloading initial positions. lam0*10 is arbitrary, the
            # aim here is to stop the system from reaching phase B.
            if (bopx > lam0*10):
                positions = initsim.initpositions(params)
                op = orderpfunc(positions, params)
                epot = totalenergyfunc(positions, params)
                fin.write('RETURNING TO PHASE A\n')
                fin.write('{0} {1} {2}\n'.format(ttot, nxtal, bopx))
                
fin.close()

# now create the dictionary with shot information.  This is pickled
# and read by shots at the subsequent interface (see takeshot.py).
shotdict = {'nshots': totalqhits,'nshotseff': totalqhits,
            'nsuccess': totalqhits,'nsuccesseff' : totalqhits,
            'successnumbers' : np.array(range(1,totalqhits+1)),
            'successweights': np.ones(totalqhits)}

# write out to pickle file
fout = open('interface0.pkl', 'wb')
pickle.dump(shotdict, fout)
fout.close()

print "I have finished up going from phase A to lambda0"
