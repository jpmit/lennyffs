# funcselector.py
# James Mithen
# j.mithen@surrey.ac.uk

""" Classes for handing out energy and mccycle functions correctly,
according to the potential.

CLASSES:
FuncSelector - handles functions for the potential energy and MC cycle
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
    NONE = 'none'
    OPTIONS = {POTENTIAL : [LEN, GAUSS],
               MCTYPE : [NVT, NPT],
               ORDERPARAM: [NTF, NLD, FRACTF, FRACLD, NONE]}

    def __init__(self, params):
        self.check_input(params)
    
    @classmethod
    def check_input(cls, params):
        """Make sure the params dictionary contains all options correctly"""
        cls.option = {}
        for oname in cls.OPTIONS:
            # check that the parameter is in the dictionary
            if oname not in params:
                raise NoInputParamError, ('{0} not found in parameter '
                                          'dictionary'.format(oname))
            # check that the option chosen is one of those allowed
            if params[oname] not in cls.OPTIONS[oname]:
                raise InputParamError, ('{0} should be one of: {1}'.\
                                        format(oname,
                                               ' '.join(cls.OPTIONS[oname])))
            # store the chosen option in the class
            cls.option[oname] = params[oname]

    @classmethod
    def TotalEnergyFunc(cls):
        """Return function that evaluates total energy"""
        if cls.option[cls.POTENTIAL] == cls.LEN:
            return energy.len_totalenergy
        elif cls.option[cls.POTENTIAL] == cls.GAUSS:
            return energy.gauss_totalenergy

    @classmethod
    def MCCycleFunc(cls):
        """Return function that computes an MC cycle"""
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
        """Return function that computes an MC cycle"""
        if cls.option[cls.ORDERPARAM] == cls.NTF:
            return orderparam.ntf
        elif cls.option[cls.ORDERPARAM] == cls.NLD:
            return orderparam.nld
        elif cls.option[cls.ORDERPARAM] == cls.FRACTF:
            return orderparam.fractf
        elif cl.option[cls.ORDERPARAM] == cls.FRACLD:
            return orderparam.fracld
        elif cl.option[cls.ORDERPARAM] == cls.NONE:
            return orderparam.default
