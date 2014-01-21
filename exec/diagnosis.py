#! /usr/bin/env python
# diagnosis.py
# 26th June 2012
# James Mithen
# 
# Analyse data at all ffs interfaces and write statistics out to file
# we write two files, 'allinterfaces.out' and 'allinterfaceseff.out'
# the latter file includes pruning, so if pruning is applied then
# the results contained in 'allinterfaceseff.out' are the correct ones
# If pruning is not applied, both files will show the same results

import readwrite
import numpy as np
from ffsfunctions import getpickparams

# get number of interfaces in simulation
params = getpickparams()
nint = params['numint']

# get flux from file times.out, which contains times to
# reach interface lambda_0
ts, op = readwrite.r2col('times.out')
flux = 1.0 / np.average(ts)

# now get all probabilities
fired = np.empty(nint)
success = np.empty(nint)
probs = np.empty(nint)
firedeff = np.empty(nint)
successeff = np.empty(nint)
probseff = np.empty(nint)

for i in range(nint):
    fin = open('interface{0}.out'.format(i+1), 'r')
    flines = fin.readlines()
    fin.close()
    # get line with fired, success numbers on
    nline = flines[2].split()
    fired[i] = float(nline[0])
    success[i] = float(nline[1])
    probs[i] = success[i] / fired[i]
    # also effectiveshots and effective success (pruning applied)
    firedeff[i] = float(nline[3])
    successeff[i] = float(nline[4])
    probseff[i] = successeff[i] / firedeff[i]

# get cumulative probabilities i.e.
# the probability we reach the final phase given that we start
# from any given interface.  NB here we are assuming that the final
# phase has an OP equal to that at the final interface
p = 1.0
cumprobs = np.zeros(nint)
cumprobseff = np.zeros(nint)
cumprobs[-1] = probs[-1]
cumprobseff[-1] = probs[-1]
for i in range(nint-1):
    cumprobs[-i - 2] = cumprobs[-i - 1] * probs[-i - 2]
    cumprobseff[-i - 2] = cumprobseff[-i - 1]*probs[-i - 2]

# critical interface is smallest interface for which P(finalphase) > 0.5
for i in range(nint):
    if cumprobs[i] > 0.5:
        critint = i
        break
for i in range(nint):
    if cumprobseff[i] > 0.5:
        critinteff = i
        break

# critical OP is OP at the critical interface
critop = params['lambdas'][critint]
critopeff = params['lambdas'][critinteff]

# get 'relative variance' of rate
# For this we use Equation (20) from:
# R.J. Allen, D. Frenkel and P ten Wolde J. Chem. Phys. 124 194111 (2006)
# a sensible estimate of the error in the rate from a single FFS simulation is
# +- rate * (relative variance)^{1/2}
relvar = np.sum((1.0 - probs)/ (probs * fired))

# write all data out to file
# see also 'allinterfaceseff.out', which includes pruning.
sname = 'allinterfaces.out'
# probability of reaching the FINAL interface, starting from lambda0
prob = np.product(probs)
# rate of forming FINAL PHASE, starting from phase A
rate = flux*prob
fout = open(sname,'w')
fout.write('SUMMARY\nCritint Critn\n{:d} {:d} {:d}\n'
           'Flux Prob Rate ln(Rate) relvar\n{:.6e} {:.6e} {:.6e} '
           '{:.6e} {:.6e}\n'.format(critint, critop, flux, prob, rate,
                                    np.log(rate), relvar))

fout.write('----------\nDETAILED BREAKDOWN\n'
           'phase A OP: {:d} lambda0 OP: {:d} phase B OP: {:d}\n'
           '{:d} shots fired from A to lambda0\n'
           'interface OP fired success P(success) P(finalphase)\n'\
           .format(params['lambdaA'], params['lambdas'][0],
                   params['lambdas'][-1], len(ts)))
fstr = ''
for i in range(nint):
    fstr = '{0}{:d}->{:d} {:d} {:d} {:d} {:.6f} {:.6e}\n'\
            .format(fstr, i, i+1, params['lambdas'][i + 1], fired[i],
                    success[i], probs[i], cumprobs[i])
fout.write(fstr)
fout.close()

# write all data out to file for pruning run
# TODO(James): remove code replication in this file
sname = 'allinterfaceseff.out'
# probability of reaching the FINAL interface, starting from lambda0
probeff = np.product(probseff)
# rate of forming FINAL phase, starting from phase A
rateeff = flux*probeff
fout = open(sname,'w')
fout.write('SUMMARY\nCritint Critn\n{:d} {:d}\nFlux Prob Rate ln(Rate)\n'
           '{:.6e} {:.6e} {:.6e} {:.6e}\n'\
           .format(critinteff, critopeff, flux, probeff, rateeff,
                   np.log(rateeff)))
fout.write('----------\nDETAILED BREAKDOWN\n'
           'phase A OP: {:.d} lambda0 OP: {:.d} phase B OP: {:.d}\n'
           '{:.d} shots fired from A to lambda0\ninterface '
           'OP fired success P(success) P(finalphase)\n'\
           .format(params['lambdaA'], params['lambdas'][0],
                   params['lambdas'][-1], len(ts)))
fstr = ''
for i in range(nint):
    fstr = '{0}{:d}->{:d} {:d} {:d} {:d} {:.6f} {:.6e}\n'\
           .format(fstr, i, i+1, params['lambdas'][i+1], firedeff[i],
                   successeff[i], probseff[i], cumprobseff[i])
fout.write(fstr)
fout.close()
