# potselector.py
# James Mithen
# j.mithen@surrey.ac.uk

""" Class for handing out energy and mccycle functions correctly,
according to the potential.

CLASSES:
PotSelector - handles functions for the potential
"""

from lenexceptions import *
import energy
import mccycle

class PotSelector(object):
    # syntactic sugar
    LEN = 'len'
    GAUSS = 'gauss'
    POTENTIALS = [LEN,  # Lennard-Jones interaction
                  GAUSS # Gaussian interaction
                  ]
    NVT = 'nvt'
    NPT = 'npt'
    CYCLETYPES = [NVT, # NVT Metropolis MC
                  NPT  # NPT Metropolis MC
                 ]
    def __init__(self, params):
        # check that params dict has keys 'potential' and 'mctype'
        if 'potential' not in params:
            raise NoInputParamError, ('potential not found in parameter'
                                      ' dictionary')
        if 'mctype' not in params:
            raise NoInputParamError, ('mctype not found in parameter'
                                      ' dictionary')
        self.potential = params['potential']
        self.cycletype = params['mctype']
        # check that potential and cycletype specified are allowed
        if self.potential not in PotSelector.POTENTIALS:
            raise InputParamError, ('potential should be one of '
                                    '%s' %(' '.join(PotSelector.POTENTIALS)))
        if self.cycletype not in PotSelector.CYCLETYPES:
            raise InputParamError, ('mctype should be one of '
                                    '%s' %(' '.join(PotSelector.POTENTIALS)))

    def TotalEnergyFunc(self):
        """Return function that evaluates total energy"""
        if self.potential == PotSelector.LEN:
            return energy.len_totalenergy
        elif self.potential == PotSelector.GAUSS:
            return energy.gauss_totalenergy
        
    def MCCycleFunc(self):
        """Return function for doing MC Cycle"""
        if self.cycletype == PotSelector.NVT:
            # functions for NVT MC for each different potential
            if self.potential == PotSelector.LEN:
                return mccycle.len_cyclenvt
            elif self.potential == PotSelector.GAUSS:
                return mccycle.gauss_cyclenvt
        if self.cycletype == PotSelector.NPT:
            # functions for NPT MC for each different potential            
            if self.potential == PotSelector.LEN:
                return mccycle.len_cyclenpt
            elif self.potential == PotSelector.GAUSS:
                return mccycle.gauss_cyclenpt
