# potselector.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Class for handing out energy and mccycle functions correctly, according to the potential
CLASSES:
PotSelector - handles functions for the potential
"""

from lenexceptions import *
import energy

class PotSelector:
    # syntactic sugar
    LEN = 'len'
    GAUSS = 'gauss'
    POTENTIALS = [LEN,  # Lennard-Jones interaction
                  GAUSS # Gaussian interaction
                  ]
    def __init__(self, params):
        if 'potential' not in params:
            raise NoInputParam, ('potential not found in parameter'
                                 'dictionary')
        self.potential = params['potential']
        if self.potential not in POTENTIALS:
            raise InputParamError, ('potential should be one of '
                                    '%s' %(' '.join(POTENTIALS)))

    def TotalEnergy(self, positions, params):
        """Return total energy of system"""
        if self.potential == LEN:
            return energy.len_totalenergy(positions, params)
        elif self.potential == GAUSS:
            return energy.gauss_totalenergy(positions, params)

    def CycleNVT(self, positions, params, etot):
        """NVT Metropolis Monte-Carlo"""
        if self.potential == LEN:
            return mccycle.len_cyclenvt(positions, params, etot)
        elif self.potential == GAUSS:
            return mccycle.gauss_cyclenvt(positions, params, etot)

    def CycleNPT(self, positions, params, etot):
        """NPT Metropolis Monte-Carlo"""
        if self.potential == LEN:
            return mccycle.len_cyclenpt(positions, params, etot)
        elif self.potential == GAUSS:
            return mccycle.gauss_cyclenpt(positions, params, etot)
