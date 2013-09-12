#! /usr/bin/env python
# code.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
This is the main Monte-Carlo (MC) code.  This code just performs a
certain number of MC cycles, governed by the parameters in the 'in'
file.  It will calculate the 'order parameter' every so often, and
save the configuration every so often (for both of these, how often is
specified in the input file. For Forward Flux Sampling (FFS)
simulations, see codeffs.py.
"""

import numpy as np
import funcselector
import initsim
import writeoutput
import energy
import mccycle
import time
import sys

class MCProgram(object):
    """The main MC program."""
    
    def __init__(self):
        """Read parameters, get initial positions, compute initial
        potential energy."""

        # read input parameters and write to file
        self.params = initsim.getparams()
        writeoutput.writepickparams(params)
        writeoutput.writeparams(params)

        # initialize positions
        self.positions = initsim.initpositions(params)

        # write initial positions to file if new simulation
        if self.params['simulation'] == 'new':
            writeoutput.writexyzld('initpositions.xyz', positions,
                                   self.params)

        # From params dictionary create FuncSelector object.  This
        # will handle correct selection of the underlying fortran/C++
        # functions correctly (the functions called depend on the
        # potential, i.e. the value of params['potential'], and also
        # on the type of MC cycle wanted, i.e.  params['mctype'], and
        # on the order parameter desired, params['orderparam'].
        funcman = funcselector.FuncSelector(params)
        self.totalenergy = funcman.TotalEnergyFunc()
        self.runcycle = funcman.MCCycleFunc()
        self.orderp = funcman.OrderParamFunc()

        # number of times to call MC cycle function
        self.ncall = int( np.ceil(params['ncycle'] /
                                  float(params['opsamp'])) )

        # number of cycles each time we call MC cycle function
        self.params['cycle'] = min(self.params['ncycle'],
                                   self.params['opsamp'])

        def run(self):
            """Perform the MC simulation."""

            # compute initial energy
            epot = totalenergy(self.positions, self.params)

            # file for writing order parameter
            opfile = open('opval.out','w')
            starttime = time.time()

            # run the MC cycles
            cyclesdone = 0
            opfile.write('{0} {1}\n'.format(0, orderp(self.positions,
                                                      self.params)))
            for cy in range(ncall):
                self.positions, epot = self.runcycle(self.positions,
                                                     self.params,
                                                     self.epot)
                cyclesdone += self.params['cycle']
                # write out order parameter
                opfile.write('{0} {1}\n'.format(cyclesdone,
                                                self.orderp(self.positions,
                                                            self.params)))
                opfile.flush()
                # write out pos file if required
                if (cyclesdone % self.params['nsave'] == 0):
                    writeoutput.writexyzld('pos{0}.xyz'.format(cy), self.positions,
                                           self.params)
        
        endtime = time.time()

        # write final positions to file
        writeoutput.writexyzld('finalpositions.xyz', self.positions,
                               self.params)

        # write runtime to stderr
        sys.stderr.write("runtime in s: %.3f\n" %(endtime - starttime))

if __name__ == '__main':
    mcprog = MCProgram()
    mcprog.run()
