#! /usr/bin/env python
# finish.py
# James Mithen
# j.mithen@surrey.ac.uk
#
# We have made it from interface intfrom to interface intfrom+1
# This script gathers the results of every shot at this interface
# At the end of the FFS run, these interfacei.out files are used
# by the script diagnosis.py, which collects FFS statistics

import sys
import os
import glob
import pickle
import numpy as np

# get arguments and complain if not right
# expect at least 1 argument: intfrom
# which is the interface we are arriving from
# e.g. if intfrom=0, we are taking a shot from lambda0, aiming for lambda1
# the second argument is optional, if set, we dont delete the files produced
# by each individual shot (this is useful for debugging, but in general we
# want to delete these files)

if len(sys.argv) != 2 and len(sys.argv) != 3:
    sys.exit('Error: finish.py expected either one or two argument')
intfrom = int(sys.argv[1]) # interface we came from

if len(sys.argv) == 3:
    delete = False
else:
    delete = True

# exit the program if we have just gone from phaseA->lambda0
# since nothing to do in this case (make this cleaner later)
# since lambda0.py writes its own shot dictionary
if (intfrom == -1) :
    sys.exit()

# work out nshots by counting number of shots%d_%d.out files
sfiles = glob.glob('shots%d_*.out' %(intfrom + 1))
nshots = len(sfiles)

# get sorted list of shot numbers
shotnums = np.zeros(nshots,dtype=int)
i = 0
for sfile in sfiles:
    snum = int(sfile.split('_')[1].split('.')[0])    
    shotnums[i] = snum
    i = i + 1
shotnums.sort()

# open each shot file for interface, and get from,time,success,weight
shotfrom = np.zeros(nshots,dtype=int)
times = np.zeros(nshots,dtype=int)
success = np.zeros(nshots,dtype=int)
weights = np.zeros(nshots)
# also want shot numbers of successful shots and their weights
successnumbers = []
successweights = []
nsuccess = 0 # total number successes
nsuccesseff = 0.0 # total weight of successful shots (>= nsuccess)
i = 0
for sn in shotnums:
    fout = open('shots%d_%d.out' %(intfrom+1,sn),'r')
    dataline = fout.readlines()[1].split()
    fout.close()
    f = int(dataline[0])
    t = int(dataline[1])
    s = int(dataline[2])
    w = float(dataline[3])
    shotfrom[i] = f
    times[i] = t
    success[i] = s
    weights[i] = w
    if s:
        successnumbers.append(sn)
        successweights.append(w)
    nsuccess = nsuccess + s
    nsuccesseff = nsuccesseff + w
    i = i + 1
# 'effective' number of shots, taking into account pruning
nshotseff = nsuccesseff + (nshots - nsuccess)

# if zero successful shots at this interface,
# then exit printing an error message
#if nsuccess == 0:
#    sys.exit('Error: no successful shots at this interface')
                           
# write shot data to interface file
# In the SUMMARY section is
# nshots - the number of shots fired to reach the interface
# nsuccess - the number of successful shots
# P(success) - fraction of shots successful
# nshotseff - the effective number of shots fired,
#             taking into account pruning (>=nfired)
# nsuccesseff - the effective number of successful shots, again
#               taking into account prunint (>=nsuccess)
# P(successeff) - fracition of effective shots successful
# Note that it is P(successeff) that is our estimate of
# reaching the next interface.
# If pruning is not applied, then P(successeff) == P(success)
# see Allen, Valerani, ten Wolde J. Phys. Condens. matter 21, 463102
# for more info on pruning.

fout = open('interface%d.out' %(intfrom+1),'w')
fout.write('SUMMARY\nshots nsuccess P(success) nshotseff nsuccesseff P(successeff)\n')
fout.write('%d %d %.6f %.6f %.6f %.6f\n'
           %(nshots,nsuccess, float(nsuccess)/float(nshots),
             nshotseff,nsuccesseff,nsuccesseff/nshotseff))
# Write a detailed breakdown of every shot
fout.write('----------\nDETAILED BREAKDOWN\nshotnum from time success\n')
fstr = ''
for (s,f,t,suc,w) in zip(shotnums,shotfrom,times,success,weights):
    # write shotnum from time success weight
    fstr = '%s%d %d %d %d %.3f\n' %(fstr,s,f,t,suc,w)
fout.write(fstr)
fout.close()

# now create the dictionary with shot information
# this is pickled and read by shots at the subsequent interface
shotdict = {'nshots': nshots,'nshotseff': nshotseff,
            'nsuccess': nsuccess,
            'nsuccesseff' : nsuccesseff,
            'successnumbers' : np.array(successnumbers),
            'successweights': np.array(successweights)}
# write out to pickle file
fout = open('interface%d.pkl' %(intfrom+1), 'wb')
pickle.dump(shotdict, fout)
fout.close()

# now clean up the directory by deleting files
# delete the error and output files produced by qsub command
# unless the script received a third argument (see above)
if delete:
    delfiles = glob.glob('shots%d*.e*' %(intfrom + 1))
    for f in delfiles:
        os.remove(f)
    delfiles = glob.glob('shots%d*.o*' %(intfrom + 1))
    for f in delfiles:
        os.remove(f)

    # delete the shots1_1.out files that are written by takeshot
    delfiles = glob.glob('shots%d*.out' %(intfrom + 1))
    for f in delfiles:
        os.remove(f)
