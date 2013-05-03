#! /usr/bin/env python
# codeffs.py
# 19th July 2012
# James Mithen
#
# FFS implementation on oracle grid engine using array job facility
# This script is an implementation of 'Forward Flux Sampling' (FFS)
# see Allen, Valerani, ten Wolde J. Phys. Condens. matter 21, 463102
# the particular algorithm implemented here is referred to as DFFS in that paper
#
# FFS is implemented in this script by submitting a lot of separate jobs

import os
import initsim
from ffsfunctions import *
import writeoutput

# path of executables lambda0.py, finish.py, takeshot.py
epath = '/user/phstf/jm0037/awork/montecarlo/epitaxy/codeffs/scripts'
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

if ffsre == 'new':
    # go from phase A to phase lambda0
    substring = ('qsub %s -cwd -N shots0_0_%s -b y '
                 '%s/lambda0.py' %(qstr,ffsnm,epath))
    os.system(substring)
    print "running command: %s" %substring
    # finish up interface
    substring = ('qsub %s -cwd -hold_jid shots0_0_%s -N finish0_%s '
                '-b y %s/finish.py -1' %(qstr,ffsnm,ffsnm,epath))
    os.system(substring)
    print "running command: %s" %substring

# for each interface in turn submit an array job
if ffsre == 'new':
    intstart = 0
else:
    intstart = int(ffsre)
    
for nint in range(intstart,numint):

    # work out job to hold for (if any)
    if nint == intstart:
        # for first interface, only hold job if we started FFS from beginning
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
                '%s/takeshot.py %d' %(qstr,holdstr, jobnm,nshots,epath,nint))
    print 'running command: %s' %substring
    os.system(substring) # old school Python way to execute string

    # take extra shots
    for bat in range(1,nbatch):
        # these extra runs just die unless we have fewer than
        # minsuccess successful shots
        jobnm = 'shots%d_%d_%s' %(nint + 1, bat, ffsnm)
        inmin = bat*nshots + 1 # min job array number
        inmax = inmin + nshots - 1 # max job array number
        substring = ('qsub %s -cwd -hold_jid shots%d_%d_%s '
                     '-N %s -t %d:%d -b y %s/takeshot.py %d yes'
                     %(qstr,nint+1,bat-1,ffsnm,jobnm,inmin,inmax,epath,nint))
        print 'running command: %s' %substring        
        os.system(substring) # old school Python way to execute string
 
    # now clean up interface
    jobnm = 'finish%d_%s' %(nint + 1, ffsnm) # pbs job name
    substring = ('qsub %s -cwd -hold_jid shots%d_%d_%s -N %s -b y ' 
                '%s/finish.py %d' %(qstr,nint + 1,nbatch-1,ffsnm,
                                    jobnm,epath,nint))
    print 'running command: %s' %substring
    os.system(substring) # old school Python way to execute string

# now that we have gone through all the interfaces, call the diagnostic code
# this will calculate FFS statistics. See file diagnosis.py.
