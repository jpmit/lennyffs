#! /usr/bin/env python
# lambda0.py
# James Mithen
# j.mithen@surrey.ac.uk
#
# Take the system from phase A to Lambda0 in FFS simulation.
# see Allen, Valerani, ten Wolde J. Phys. Condens. matter 21, 463102

import sys
import pickle
import numpy as np
import initsim
import energy
import mccycle
from bops import bopxbulk
from ffsfunctions import savelambda0config, getpickparams
from writeoutput import  writexyz
import potselector
from lenexceptions import *

params = getpickparams()

lamA = params['lambdaA']
lam0 = params['lambdas'][0]
totalqhits = params['totalqhits'] # num times to go through lambda0
                                  # from phase A
lamsamp = params['lambdasamp'] # number of MC cycles per OP evaluation
params['cycle'] = lamsamp

# initialize positions
positions = initsim.initpositions(params)

# check that lambda < lamA (we are in phase A)
nxtal, bopx = bopxbulk(positions,params)
if (bopx >= lamA):
    raise FFSError, ('BOPxbulk is {0}, system must start in phase A, '
                     'BOPxbulk < {1}'.format(bopx,lamA))

# get the correct energy function and MC cycle function using
# PotSelector interface

PotManager = potselector.PotSelector(params)
totalenergyfunc = PotManager.TotalEnergyFunc()
cyclefunc = PotManager.MCCycleFunc()

# write initial OP to file
fin = open('ffsacc.out','w')
fin.write('ncycles nxtal bopX\n')
fin.write('%d %d %d\n' %(0,nxtal,bopx))
fin.flush()

# start simulation
qhits = 0 # num times passed through lambda0
thit = 0
ttot = 0

epot = totalenergyfunc(positions,params)    
while (qhits < totalqhits):
    # make some trial moves
    positions,epot = cyclefunc(positions,params,epot)
    thit = thit + lamsamp
    ttot = ttot + lamsamp
    # evaluate OP
    nxtal,bopx = bopxbulk(positions,params)
    # write OP and total time to file
    fin.write('%d %d %d\n' %(ttot,nxtal,bopx))
    fin.flush()
    if (bopx >= lam0):
        # we hit the interface
        qhits = qhits + 1
        savelambda0config(qhits,thit,positions,params)
        thit = 0
        # now let the system relax back to lambda A
        while (bopx >= lamA):
            positions,epot = cyclefunc(positions,params,epot)
            ttot = ttot + lamsamp                
            # evaluate OP
            nxtal, bopx = bopxbulk(positions,params)
            # write OP and total time to file
            fin.write('%d %d %d\n' %(ttot,nxtal,bopx))
            # if we reached lam0*10 (arbitrary, means going to phase B),
            # then return to phase A
            if (bopx > lam0*10):
                positions = initsim.initpositions(params)
                nxtal,bopx = bopxbulk(positions,params)
                epot = totalenergyfunc(positions,params)
                fin.write('RETURNING TO PHASE A\n')
                fin.write('%d %d %d\n' %(ttot,nxtal,bopx))
                
fin.close()

# now create the dictionary with shot information
# this is pickled and read by shots at the subsequent interface
# (see takeshot.py)
shotdict = {'nshots': totalqhits,'nshotseff': totalqhits,
            'nsuccess': totalqhits,'nsuccesseff' : totalqhits,
            'successnumbers' : np.array(range(1,totalqhits+1)),
            'successweights': np.ones(totalqhits)}
# write out to pickle file
fout = open('interface0.pkl', 'wb')
pickle.dump(shotdict, fout)
fout.close()

print "I have finished up going from phase A to lambda0"
