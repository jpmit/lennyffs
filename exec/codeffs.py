#! /usr/bin/env python
# codeffs.py
# 26th June 2012
# James Mithen
#
# FFS implementation on oracle grid engine using array job facility
# This script is an implementation of 'Forward Flux Sampling' (FFS)
# see Allen, Valerani, ten Wolde J. Phys. Condens. matter 21, 463102.
# The particular algorithm implemented here is referred to as DFFS in
# that paper.
#
# FFS is implemented in this script by submitting a lot of separate
# jobs to the cluster, using the 'qsub' command.  Each job is in fact
# an 'array' job, meaning that each job is essentially a wrapper for
# a large number of jobs.  The number of jobs in the array is the
# number of shots fired at the interface.
#
# For further information
# about the qsub command, array jobs, and
# other options, see the oracle grid engine manual.  It is hoped
# that this script could be modified to be run with other job
# submission systems.
#
# Jobs for all of the interfaces are submitted at the same time:
# the 'hold' facility is used so that jobs at e.g. interface 2 won't
# start running until the jobs at interface 1 have all completed.
#
# Also, there is a facility to run a variable number of shots at each
# interface.  'nshots', specified in the input file, is the number
# of shots initially run for each interface.  All of these shots are
# guaranteed to be taken.  If these shots result in at least
# 'minsuccess' successes (i.e. shots that reach the next interface),
# no more shots are taken.  If not, additional batches of shots are
# taken in turn, until the minimum number of successes are reached.
# Note that only 'nbatch' batches of shots will be taken in total,
# so the maximum number of shots that will be taken is nbatch*nshots.
#
# After all of the batches of shots have finished for a particular
# interface, the script 'finish.py' is run.  This script will clean
# up the working directory, and output a file 'interface[inum].out',
# which has details of the number of shots taken, the ones that were
# successfull, etc.
#
# Once all the jobs submitted by this script have finished, you need
# to run 'diagnosis.py' in this directory.  This will read the relevant
# data from the interface*.out files and compute FFS statistics for the
# entire simulation i.e. rates.

import os
import initsim
from ffsfunctions import *
import writeoutput

# path of executables lambda0.py, finish.py, takeshot.py
# this is currently set to be ./../scripts relative to this directory
# there is a bit of a hack here to step back a directory
epath = os.path.expanduser('~/awork/montecarlo/epitaxy/codeffs/scripts/')

# can specify a job queue here
qname = ''
if qname:
    qstr = '-q %s' %qname
else:
    qstr = ''

# get params and write to pickle file 'params.pkl' for future reading
params = initsim.getparams()
writeoutput.writepickparams(params)

# write params to 'pickle.out' -> human readable version of params.pkl
writeoutput.writeparams(params)

# params needed for job submission
numint = params['numint']
nshots = params['nshots']
nbatch = params['nbatch']
minsuccess = params['minsuccess']
ffsnm = params['ffsname']
ffsre = params['ffsrestart']

# if ffsrestart in params file is set to -1, we start a 'new' FFS
# simulation
if ffsre == -1:
    # go from phase A to phase lambda0
    substring = ('qsub %s -cwd -N shots0_0_%s -b y '
                 '%s/lambda0.py' %(qstr,ffsnm,epath))
    os.system(substring) # old school Python way to execute string
    print "running command: %s" %substring
    # finish up interface
    substring = ('qsub %s -cwd -hold_jid shots0_0_%s -N finish0_%s '
                '-b y %s/finish.py -1' %(qstr,ffsnm,ffsnm,epath))
    os.system(substring)
    print "running command: %s" %substring

# for each interface in turn submit an array job
if ffsre == -1:
    intstart = 0
else:
    intstart = ffsre
    
for nint in range(intstart,numint):

    # work out job to hold for (if any)
    if nint == intstart:
        # for first interface, only hold job if we started
        #  FFS from beginning
        if ffsre == 'new':
            holdstr = '-hold_jid finish0_%s' %ffsnm
        else:
            holdstr = ''
    else:
        # hold for finish of previous interface
        holdstr = '-hold_jid finish%d_%s' %(nint,ffsnm)

    # take shots
    jobnm = 'shots%d_0_%s' %(nint + 1,ffsnm) # pbs job name        
        
    substring = ('qsub %s -cwd %s -N %s -t 1:%d -b y ' 
                '%s/takeshot.py %d' %(qstr,holdstr, jobnm,nshots,
                                      epath,nint))
    print 'running command: %s' %substring
    os.system(substring)

    # take extra shots
    for bat in range(1,nbatch):
        # these extra runs just die unless we have fewer than
        # minsuccess successful shots
        jobnm = 'shots%d_%d_%s' %(nint + 1, bat, ffsnm)
        inmin = bat*nshots + 1 # min job array number
        inmax = inmin + nshots - 1 # max job array number
        substring = ('qsub %s -cwd -hold_jid shots%d_%d_%s '
                     '-N %s -t %d:%d -b y %s/takeshot.py %d yes'
                     %(qstr,nint+1,bat-1,ffsnm,jobnm,inmin,inmax,
                       epath,nint))
        print 'running command: %s' %substring        
        os.system(substring)
 
    # now clean up interface via finish.py script
    jobnm = 'finish%d_%s' %(nint + 1, ffsnm) # pbs job name
    substring = ('qsub %s -cwd -hold_jid shots%d_%d_%s -N %s -b y ' 
                '%s/finish.py %d' %(qstr,nint + 1,nbatch-1,ffsnm,
                                    jobnm,epath,nint))
    print 'running command: %s' %substring
    os.system(substring) 
