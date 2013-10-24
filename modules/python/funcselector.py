# funcselector.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Class for handing out energy and mccycle functions correctly,
according to the potential.

CLASSES:
FuncSelector - handles function selection for the potential energy,
               MC cycle and order parameter.
"""

from lenexceptions import *
import energy
import mccycle
import orderparam

class FuncSelector(object):
    """Interface for selecting correct functions for:
    - Interaction Potential
    - MC Cycle type
    - Order parameter
    """
    POTENTIAL = 'potential'
    MCTYPE = 'mctype'
    # choices for potential
    LEN = 'len'
    GAUSS = 'gauss'
    # choices for mctype
    NVT = 'nvt'
    NPT = 'npt'
    # choices for orderparam
    ORDERPARAM = 'orderparam'
    NTF = 'ntf'
    NLD = 'nld'
    FRACTF = 'fractf'
    FRACLD = 'fracld'
    ALLFRACLD = 'allfracld'
    ALLFRAC = 'allfrac'
    NONE = 'none'
    # note the first option is taken as the default in each case
    OPTIONS = {POTENTIAL : [LEN, GAUSS],
               MCTYPE : [NVT, NPT],
               ORDERPARAM: [NTF, NLD, FRACTF, FRACLD, ALLFRACLD,
                            ALLFRAC, NONE]}

    def __init__(self, params):
        self.store_input(params)
    
    @classmethod
    def store_input(cls, params):
        """
        Store the chosen options for interaction potential, MC cycle
        type, and order parameter in the instance.
        """
        cls.option = {}
        for oname in cls.OPTIONS:
            # store the chosen option in the class
            try:
                cls.option[oname] = params[oname]
            except:
                # didn't get option: default to the first one
                cls.option[oname] = cls.OPTIONS[oname][0]
                pass

    @classmethod
    def TotalEnergyFunc(cls):
        """Return function that evaluates total energy."""
        if cls.option[cls.POTENTIAL] == cls.LEN:
            return energy.len_totalenergy
        elif cls.option[cls.POTENTIAL] == cls.GAUSS:
            # this uses neighbour lists
            return energy.gauss_totalenlist

    @classmethod
    def EnergyIparFunc(cls):
        """Return function that evaluates energy of a single particle."""
        if cls.option[cls.POTENTIAL] == cls.LEN:
            return energy.len_energyipar
        elif cls.option[cls.POTENTIAL] == cls.GAUSS:
            return energy.gauss_energyipar

    @classmethod
    def MCCycleFunc(cls):
        """Return function that computes an MC cycle."""
        if cls.option[cls.MCTYPE] == cls.NPT:
            # functions for NPT MC for each different potential            
            if cls.option[cls.POTENTIAL] == cls.LEN:
                return mccycle.len_cyclenpt
            elif cls.option[cls.POTENTIAL] == cls.GAUSS:
                return mccycle.gauss_cyclenpt
        if cls.option[cls.MCTYPE] == cls.NVT:
            # functions for NVT MC for each different potential            
            if cls.option[cls.POTENTIAL] == cls.LEN:
                return mccycle.len_cyclenvt
            elif cls.option[cls.POTENTIAL] == cls.GAUSS:
                return mccycle.gauss_cyclenvt

    @classmethod
    def OrderParamFunc(cls):
        """Return function that computes an MC cycle."""
        if cls.option[cls.ORDERPARAM] == cls.NTF:
            return orderparam.nclustf_cpp
        elif cls.option[cls.ORDERPARAM] == cls.NLD:
            return orderparam.nclusld_cpp
        elif cls.option[cls.ORDERPARAM] == cls.FRACTF:
            return orderparam.fractf_cpp
        elif cls.option[cls.ORDERPARAM] == cls.FRACLD:
            return orderparam.fracld_cpp
        elif cls.option[cls.ORDERPARAM] == cls.ALLFRACLD:
            return orderparam.allfracld_cpp
        elif cls.option[cls.ORDERPARAM] == cls.ALLFRAC:
            return orderparam.allfracldtf_cpp        
        elif cls.option[cls.ORDERPARAM] == cls.NONE:
            return orderparam.default
