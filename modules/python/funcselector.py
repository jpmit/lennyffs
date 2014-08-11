# funcselector.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Class for handing out energy and mccycle functions correctly,
according to the potential.  At the moment the code the Lennard-Jones
and Gaussian potentials only.

CLASSES:
FuncSelector - handles function selection for the potential energy,
               MC cycle and order parameter.
"""

from lenexceptions import *
import energy
import mccycle
import orderparam
import writeoutput

class FuncSelector(object):
    """
    Interface for selecting correct functions for:
    - Interaction Potential
    - MC Cycle type
    - Order parameter
    """
    
    POTENTIAL = 'potential'
    MCTYPE = 'mctype'
    ORDERPARAM = 'orderparam'
    WRITEXYZ = 'writexyz'
    # choices for potential
    LEN = 'len'
    GAUSS = 'gauss'
    IPL = 'ipl'
    # choices for mctype
    NVT = 'nvt'
    NPT = 'npt'
    MD  = 'md'
    # choices for orderparam
    Q6 = 'q6global' # global Q6 of entire system
    NTF = 'ntf' # largest cluster according to TF
    NLD = 'nld' # largest cluster according to LD
    FRACTF = 'fractf' # frac xtal particles according to TF
    FRACLD = 'fracld' # frac xtal particles according to LD
    ALLFRACLD = 'allfracld' # frac of all polymorphs (LD)
    ALLFRAC = 'allfrac' # frac of all polymorphs (LD) + TF xtal frac
    POTENERGY = 'potenergy' # potential energy of the system
    NONE = 'none'
    # choices for writexyz
    LD = 'ld'
    TF = 'tf'
    NOOP = 'noop'
    # the first on the list is taken as the default here (!)
    OPTIONS = {POTENTIAL : [LEN, GAUSS],
               MCTYPE : [NVT, NPT, MD],
               ORDERPARAM: [Q6, NTF, NLD, FRACTF, FRACLD, ALLFRACLD,
                            ALLFRAC, NONE],
               WRITEXYZ: [TF, LD, NOOP]
               }

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
            # this uses neighbour lists
            return energy.len_totalenlist
        elif cls.option[cls.POTENTIAL] == cls.GAUSS:
            # this uses neighbour lists
            return energy.gauss_totalenlist
        elif cls.option[cls.POTENTIAL] == cls.IPL:
            return energy.ipl_totalenlist

    @classmethod
    def EnergyIparFunc(cls):
        """Return function that evaluates energy of a single particle."""
        
        if cls.option[cls.POTENTIAL] == cls.LEN:
            return energy.len_energyipar
        elif cls.option[cls.POTENTIAL] == cls.GAUSS:
            return energy.gauss_energyipar
        elif cls.option[cls.POTENTIAL] == cls.IPL:
            return energy.ipl_energyipar

    @classmethod
    def MCCycleFunc(cls):
        """Return function that computes an MC cycle."""
        
        if cls.option[cls.MCTYPE] == cls.NPT:
            # functions for NPT MC for each different potential            
            if cls.option[cls.POTENTIAL] == cls.LEN:
                return mccycle.len_cyclenpt
            elif cls.option[cls.POTENTIAL] == cls.GAUSS:
                return mccycle.gauss_cyclenpt
            elif cls.option[cls.POTENTIAL] == cls.IPL:
                return mccycle.ipl_cyclenpt
        if cls.option[cls.MCTYPE] == cls.NVT:
            # functions for NVT MC for each different potential            
            if cls.option[cls.POTENTIAL] == cls.LEN:
                return mccycle.len_cyclenvt
            elif cls.option[cls.POTENTIAL] == cls.GAUSS:
                return mccycle.gauss_cyclenvt
            elif cls.option[cls.POTENTIAL] == cls.IPL:
                return mccycle.ipl_cyclenvt
        if cls.option[cls.MCTYPE] == cls.MD:
            # NVE MD is only available for Gaussian potential currently
            if cls.option[cls.POTENTIAL] == cls.GAUSS:
                return mccycle.gauss_cyclemd

    @classmethod
    def OrderParamFunc(cls):
        """Return function that computes an MC cycle."""

        if cls.option[cls.ORDERPARAM] == cls.Q6:
            return orderparam.q6global_cpp
        elif cls.option[cls.ORDERPARAM] == cls.NTF:
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
        elif cls.option[cls.ORDERPARAM] == cls.POTENERGY:
            # potential energy OP depends on the interaction type
            if cls.option[cls.POTENTIAL] == cls.LEN:
                # output potential energy per particle in units of
                # epsilon, rather than in units of 4 epsilon
                def len_totalen_epsilon(positions, params):
                    return 4.0 * energy.len_totalenlist(positions, params) / params['npartot']
                return len_totalen_epsilon
            elif cls.option[cls.POTENTIAL] == cls.GAUSS:
                def gauss_totalen_epsilon(positions, params):
                    return energy.gauss_totalenlist(positions, params) / params['npartot']
                return gauss_totalen_epsilon
            elif cls.option[cls.POTENTIAL] == cls.IPL:
                def ipl_totalen_epsilon(positions, params):
                    return energy.ipl_totalenlist(positions, params) / params['npartot']
                return ipl_totalen_epsilon
        elif cls.option[cls.ORDERPARAM] == cls.NONE:
            return orderparam.nothing

    @classmethod
    def WriteXyzFunc(cls):
        """Return function that will write an XYZ file."""
        
        if cls.option[cls.WRITEXYZ] == cls.LD:
            return writeoutput.writexyz_ld
        elif cls.option[cls.WRITEXYZ] == cls.TF:
            return writeoutput.writexyz_tf
        elif cls.option[cls.WRITEXYZ] == cls.NOOP:
            return writeoutput.writexyz_noop
