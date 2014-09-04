# James Mithen

"""Code for converting umbrella sampling histograms into a free energy curve for the 2-d case."""

import glob
import sys
from itertools import izip

import numpy as np
import scipy.optimize

import readwrite

def free_energy(nvals, aks, b, order):
    """Free energy for a single window (2d extension of Auer and Frenkel)."""

    n1s = np.array([n[0] for n in nvals])
    n2s = np.array([n[1] for n in nvals])

    result = np.zeros(len(n1s))
    index = 0
    for i in range(order + 1):
        n1s_power = n1s**i
        for j in range(order + 1 - i):
            #result += aks[index] * n1s**i * n2s**j
            result += aks[index] * n1s_power * n2s**j
            index += 1
    
    return result + b

def free_energy_multiple(nvals, aks, bis, order):

    # note each window can be a different length (have a different number of points)
    prediction = np.array([])
    for w, nval in enumerate(nvals): # w is window number
        prediction = np.append(prediction, free_energy(nvals[w], aks, bis[w], order))

    return prediction


def f_multiple(p, nvals, gvals, order):
    
    # number of windows nvals should be a list
    numw = len(nvals)

    bis = p[:numw]
    aks = p[numw:]

    return (gvals - free_energy_multiple(nvals, aks, bis, order))


class Stitcher2(object):
    def __init__(self):
        pass

    def load_histfiles(self, histfiles):
        """Load a set of histogram files."""

        nvals = [] # each entry gives nvals for
        gvals = []
        allgvals = np.array([]) # 1-d array

        # each list in self._nvals is values [(n1, n2), (n1, n2),...]
        # for a single window
        self._nvals = []
        self._gvals = []

        self._allgvals = np.array([])
        #print histfiles
        for i, f in enumerate(histfiles):
            n1s, n2s, gvs = readwrite.r3col(f)
            #print n1s
            ns = zip(n1s, n2s)
            self._nvals.append(ns)
            self._gvals.append(gvs)
            self._allgvals = np.append(self._allgvals, gvs)
        
        return ns, gvs

    def load_files(self, folders):
        """Load histogram data from all fhist*.out files in each of the folders."""

        # clear all existing data
        self._data = {}
        for folder in folders:
            print folder
            histfiles = glob.glob('{}/fhist*'.format(folder))

            # store the data
            self._data[folder] = self.load_histfiles(histfiles)

        # store the number of windows
        self._nw = len(histfiles)

    def fit(self, order=2):

        # order of polynomial fit
        self._order = order
        
        bistart = [0.0] * self._nw

        # for polynomial fit of order k we have ((k + 1)*(k + 2) / 2) terms
        akstart = [0.0] * ((order + 1)*(order + 2) / 2)

        # initial estimate of the parameters vector
        p0 = bistart + akstart

        p, success = scipy.optimize.leastsq(f_multiple, p0, args=(self._nvals, self._allgvals, self._order))
        self._params = p

    def get_curve(self):
        bivals = self._params[:self._nw]
        akvals = self._params[self._nw:]
        
        fecurve = np.array([])
        alln1vals = np.array([])
        alln2vals = np.array([])
        for i, nv in enumerate(self._nvals):
            bi = bivals[i]
            fe = free_energy(nv, akvals, bi, self._order)
            fe -= bi
            fecurve = np.append(fecurve, fe)

            alln1vals = np.append(alln1vals, [n[0] for n in nv])
            alln2vals = np.append(alln2vals, [n[1] for n in nv])            

        # note the shift applied to the free energy curve here
        fecurve = fecurve - fecurve[0]

        # sort since some of the windows will overlap
        sorted_lists = sorted(izip(alln1vals, alln2vals, fecurve), key = lambda x: x[0])
        alln1vals, alln2vals, fecurve = [[s[i] for s in sorted_lists] for i in range(3)]

        return alln1vals, alln2vals, fecurve


def write_output(fname, xs, ys, zs):
    outf = open(fname, 'w')
    val = xs[0]
    fstr = ''
    for (x, y, z) in zip(xs, ys, zs):
        if (x != val):
            val = x
            fstr = '{}\n'.format(fstr)
        fstr = '{}{} {} {}\n'.format(fstr, x, y, z)
    outf.write(fstr)
    outf.close()
    

if __name__ == "__main__":
    import time
    import os
    from sys import argv
    from os import getcwd
    starttime = time.time()

    # to use this program from the command line do:
    # python stitch.py [keyword1] [arg1] [keyword2] [arg2]...
    # current keywords are: folders --> ARG: 'a glob here' (must have quotes)
    #                       order   --> ARG: integer (the order of the fit)


    if 'folders' in argv:
        folderglob = argv[argv.index('folders')+1]
    else:
        folderglob = os.getcwd()

    if 'order' in argv:
        fit_order = int(argv[argv.index('order')+1])
    else:
        fit_order = 3

    folders  = glob.glob(folderglob)
             
    stitch = Stitcher2()

    r = stitch.load_files(folders)

    # order=3 gives all terms up to and including 
    stitch.fit(order=fit_order)

    n1, n2, fe = stitch.get_curve()

    # write data to output file for plotting with gnuplot
    write_output('dataout.out', n1, n2, fe)

    endtime = time.time()
    print "Total run time: " + str(endtime-starttime) + 's'
