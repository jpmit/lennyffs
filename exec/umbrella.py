#! /usr/bin/env python
# code.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
This is the umbrella-enabled main Monte-Carlo (MC) code.
This code just performs a certain number of MC cycles, governed by the
parameters in the 'in' file.  It will calculate the 'order parameter'
every so often, and save the configuration every so often (for both of these,
how often is specified in the input file). Makes use of umbrella sampling

Run with -w [index] for optional indexed output files. Index can be
anything
"""

import numpy as np
import funcselector
import initsim
import writeoutput
import energy
import force
import mccycle
import time
import sys
from copy import deepcopy
import random
import os
import pickle

class MProgram(object):
    """The main MC/MD program."""
    
    def __init__(self):
        """
        Read parameters, get initial positions,
        and setup any parameters needed for the simulation.
        """
        #allow taking index as command line argument
        if len(sys.argv) != 1 and sys.argv[1] == '-w':
            self.iwind = sys.argv[2]
        else:
            self.iwind = ''
        # read input parameters and write to file
        if self.iwind != '' and os.path.exists('params{0}.pkl'.format(self.iwind)):
            self.params = pickle.load(open('params{0}.pkl'.format(self.iwind)))
        else:
            self.params = initsim.getparams()
        if self.iwind == '':  #indexed pkl files wil already exist otherwise
            # pickled version 'params.pkl'
            writeoutput.writepickparams(self.params)
            # human readable version 'params.out'
            writeoutput.writeparams(self.params)

        # From params dictionary create FuncSelector object.  This
        # will handle correct selection of the underlying fortran/C++
        # functions correctly (the functions called depend on the
        # potential, i.e. the value of params['potential'], and also
        # on the type of MC cycle wanted, i.e.  params['mctype'], and
        # on the order parameter desired, params['orderparam'].
        funcman = funcselector.FuncSelector(self.params)
        self.totalenergy = funcman.TotalEnergyFunc()
        self.runcycle = funcman.MCCycleFunc()
        self.orderp = funcman.OrderParamFunc()
        self.writexyz = funcman.WriteXyzFunc()

        # initialize positions
        self.positions = initsim.initpositions(self.params)

        # write initial positions to file if new simulation
        if self.params['simulation'] == 'new':
            self.writexyz('initpositions{0}.xyz'.format(self.iwind), self.positions,
                          self.params)

        # number of times to call MC cycle function
        self.numbrellacycles = self.params['numbrellacycles']
         

        # number of cycles each time we call MC cycle function
        self.params['cycle'] = self.params['nunbiased']

    def run(self):
        """Perform the MC simulation."""
        self.run_mc()

    def run_mc(self):
        """Perform the MC simulation."""

        # compute initial energy
        epot = self.totalenergy(self.positions, self.params)

        # file for writing order parameter
        opfile = open('opval{0}.out'.format(self.iwind),'w')

        # write time 0 and initial orderp
        cyclesdone = 0
        opfile.write('{0} {1}\n'.format(0, self.orderp(self.positions,
                                                       self.params)))
        #Initialise w and N
        self.N0 = self.params['N0']
        self.N = self.orderp(self.positions,self.params)
        self.w = wfunc(self.N,self.N0,self.params['k'])
        
        starttime = time.time()

        for cy in range(self.numbrellacycles):

            # Store values that may be reverted if bias-chain is rejected
            temppositions = deepcopy(self.positions)
            templboxx, templboxy, templboxz = self.params['lboxx'], self.params['lboxy'], self.params['lboxz']
            tempepot = epot
            tempw = self.w
            tempN = self.N
            
            self.positions, epot = self.runcycle(self.positions,
                                                 self.params,
                                                 epot)

            # w test and revert to temp values if rejected
            self.N = self.orderp(self.positions,self.params)
            self.w = wfunc(self.N,self.N0,self.params['k'])
            biasprob = min(1.0,np.exp(-1.0*(self.w-tempw)))
            random.seed()
            temprand = random.random()

            # uncomment next 5 lines for testing
            #print '--------------------------------------------------'
            #print 'self.N',self.N,'tempN',tempN
            #print 'self.w',self.w,'tempw',tempw
            #print 'temprand',temprand,'biasprob',biasprob
            #print '--------------------------------------------------'
            
            if temprand > biasprob:
                self.positions = deepcopy(temppositions)
                self.params['lboxx'], self.params['lboxy'], self.params['lboxz'] = templboxx, templboxy, templboxz
                self.w = tempw
                self.N = tempN
                epot = tempepot
            
            cyclesdone += self.params['cycle']
            # write out order parameter
            opfile.write('{0} {1}\n'.format(cyclesdone,
                                            self.orderp(self.positions,
                                                        self.params)))
            opfile.flush()
            # write out pos file if required
            if (cyclesdone % self.params['nsave'] == 0):
                self.writexyz('{0}pos{1}.xyz'.format(self.iwind,cyclesdone), self.positions,
                              self.params)

        endtime = time.time()

        # write final positions to file
        self.writexyz('finalpositions{0}.xyz'.format(self.iwind), self.positions, self.params)

        # write runtime to stderr
        sys.stderr.write("runtime in s: {:.3f}\n".format(endtime - starttime))

        # if we were npt, print new box volume
        if self.params['mctype'] == 'npt':
            sys.stderr.write('new box volume: {0}\n'\
                             .format(self.params['lboxx']*\
                                     self.params['lboxy']*\
                                     self.params['lboxz']))

#------------------------------------------------------
#Bias function definition
def wfunc(N,N0,k):
    return  0.5*k*(N-N0)**2
#------------------------------------------------------

if __name__ == '__main__':
    mcprog = MProgram()
    mcprog.run()
