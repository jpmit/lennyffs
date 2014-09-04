# stitch.py
# James Mithen

"""Code for converting umbrella sampling histograms into a free energy curve for the 1-d case."""

import sys
from itertools import izip
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

def free_energy(nvals, aks, b):
    """Free energy for a single window (see Auer and Frenkel)."""

    k = len(aks) # order of polynomial fit
    
    return (np.sum([aks[i] * nvals**(i + 1) for i in range(k)], axis=0) + b)

def free_energy_multiple(nvals, aks, bis):

    # note each window can be a different length (have a different number of points)
    prediction = np.array([])
    for w, nval in enumerate(nvals): # w is window number
        prediction = np.append(prediction, free_energy(nvals[w], aks, bis[w]))

    return prediction


def f_multiple(p, nvals, gvals, weights):
    
    # number of windows nvals should be a list
    numw = len(nvals)

    bis = p[:numw]
    aks = p[numw:]

    return weights**0.5 * (gvals - free_energy_multiple(nvals, aks, bis))


class Stitcher(object):
    def __init__(self):
        pass

    def load_histfiles(self, histfiles):
        """Load a set of histogram files."""

        nvals = [] # each entry gives nvals for
        gvals = []
        allnvals = np.array([]) # 1-d array
        allgvals = np.array([]) # 1-d array

        for i, f in enumerate(histfiles):
            nvs, gvs = readwrite.r2col(f)
            if (i == 0):
                # chop of everything before the minimum for first window
                minindx = list(gvs).index(min(gvs))
            else:
                minindx = 0
            nvals.append(nvs[minindx:])
            gvals.append(gvs[minindx:])            
            allnvals = np.append(allnvals, nvals[minindx:])
            allgvals = np.append(gvals, gvs[minindx:])

        return nvals, gvals, allnvals, allgvals

    def load_files(self, folders):
        """Load histogram data from all fhist*.out files in each of the folders."""

        # clear all existing data
        self._data = {}

        sorter = lambda x : int(''.join([i for i in x if i.isdigit()]))

        nw = None # number of windows
        for fol in folders:
            histfiles = glob.glob('{}/fhist*'.format(fol))
            # sort these files
            histfiles.sort(key=sorter)

            newnw = len(histfiles)
            if (nw is not None) and (newnw != nw):
                sys.exit('incorrect number of windows for folder {}'.format(fol))
            else:
                nw = newnw

            # store the data
            self._data[fol] = self.load_histfiles(histfiles)

        # store the number of windows
        self._nw = nw

        # compute nvals, gvals, and variance in gvals
        self.compute_windows()

    def compute_windows(self):
        # get 'averaged' windows (averaged over all folders)
        self._nvals = []
        self._gvals = []
        self._gvar = []
        self._weights = []
        for w in range(self._nw):
            gvs = {}
            for k in self._data:
                # values of n and g for window w and this specific folder
                ns = self._data[k][0][w]
                gs = self._data[k][1][w]                
                for (n, g) in zip(ns, gs):
                    if n not in gvs:
                        gvs[n] = []
                    gvs[n].append(g)
            # now get the n values, g values, and variances for this window
            nvals = gvs.keys()
            nvals.sort()
            self._nvals.append(nvals)
            self._gvals.append([])
            self._gvar.append([])
            self._weights.append([])            
            for n in nvals:
                gvals = gvs[n]
                self._gvals[-1].append(np.average(gvals))
                self._gvar[-1].append(np.var(gvals))
            # go through each measurement of the variance; if the
            # variance is zero, this probably means that we only have
            # a single measurement of the fe for this cluster size.
            # In this case, we set the variance to be equal to the
            # maximum variance in the window.  If we only have a
            # single measurement, we set the variance to be 1.0 for
            # each window, which will give us weights of 1.0 in the
            # GLS fit.
            maxvar = max(self._gvar[-1])
            for i in range(len(self._gvar[-1])):
                if self._gvar[-1][i] == 0.0:
                    if maxvar == 0.0: # must have had a single folder
                        self._gvar[-1][i] = 1
                    else:
                        self._gvar[-1][i] = maxvar
            # set weights for window as inverse of variance
            self._weights[-1] = 1.0 / np.array(self._gvar[-1])
            self._nvals[-1] = np.array(self._nvals[-1])
            self._gvals[-1] = np.array(self._gvals[-1])
        self._allweights = np.array([])
        self._allgvals = np.array([])
        self._allnvals = np.array([])        
        for w in range(self._nw):
            self._allweights = np.append(self._allweights, self._weights[w])
            self._allgvals = np.append(self._allgvals, self._gvals[w])
            self._allnvals = np.append(self._allnvals, self._nvals[w])            

    def fit(self, useweights=True, order=4):
        
        bistart = [0.0] * self._nw
        akstart = [0.0] * order # number of entries is order of polynomial

        # initial estimate of the parameters vector
        p0 = bistart + akstart

        if useweights:
            weights = self._allweights
            #print weights
        else:
            weights = np.ones(len(self._allgvals))

        p, success = scipy.optimize.leastsq(f_multiple, p0, args=(self._nvals, self._allgvals, weights))
        #print success
        self._params = p

    def get_curve(self):
        bivals = self._params[:self._nw]
        akvals = self._params[self._nw:]
        fecurve = np.array([])
        allnvals = np.array([])
        for i, nv in enumerate(self._nvals):
            bi = bivals[i]
            fe = free_energy(nv, akvals, bi)
            fe -= bi
            fecurve = np.append(fecurve, fe)

            allnvals = np.append(allnvals, nv)

        # note the shift applied to the free energy curve here
        fecurve = fecurve - fecurve[0]

        # sort since some of the windows will overlap
        sorted_lists = sorted(izip(allnvals, fecurve), key=lambda x: x[0])
        ns, gs = [[x[i] for x in sorted_lists] for i in range(2)]
        return np.array(ns), np.array(gs)
        #return allnvals, fecurve - fecurve[0]


if __name__ == "__main__":

    import readwrite
    import glob
    from sys import argv
    from os import getcwd

    # to use this program from the command line do:
    # python stitch.py [keyword1] [arg1] [keyword2] [arg2]...
    # current keywords are: folders --> ARG: 'a glob here' (must have quotes)
    #                       order   --> ARG: integer (the order of the fit)
    #                       mode    --> ARG: 'separate' or 'combined' for how the
    #                                        curves are fit

    if 'folders' in argv:
        folderglob = argv[argv.index('folders')+1]
    else:
        folderglob = os.getcwd()

    if 'order' in argv:
        fit_order = int(argv[argv.index('order')+1])
    else:
        fit_order = 4

    if 'mode' in argv:
        mode = argv[argv.index('mode')+1]
    else:
        mode = 'separate'

    folders = glob.glob(folderglob)
    
    stitch = Stitcher()

    if mode == 'combined':
        stitch.load_files(folders[:2])
        stitch.fit(useweights=True)
        ns, gs = stitch.get_curve()
        plt.plot(ns, gs, label='combined with weights', color='k', linestyle='--')
        stitch.fit(useweights=False, order=fit_order)
        ns, gs = stitch.get_curve()
        plt.plot(ns, gs, label='combined', color='k')

    if mode == 'separate':
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '0.25', '0.75', '0.5']

        maxs = []
        for i, f in enumerate(folders):
            stitch.load_files([f])
            stitch.fit(useweights=True, order=fit_order)
            ns, gs = stitch.get_curve()
            print np.amax(gs)
            maxs.append(np.amax(gs))        
            #handns, handgs = readwrite.r2col(f + '/finalstitch.out')
            plt.plot(ns, gs, color=colors[i], label=f.split('/')[-1])
            #plt.plot(handns, handgs, color=colors[i], linestyle='--')  
        print 'Ave: ' + str(np.mean(maxs)) + ' +/- ' + str(np.std(maxs)) 

    plt.legend()

    plt.xlabel('$n_{ld}$')
    plt.ylabel('$G / k_b T$')
    plt.show()
