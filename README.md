LennyFFS
========

Copyright (c) 2012-2014 James Mithen.
j.mithen@surrey.ac.uk

INTRODUCTION
------------

LennyFFS is a code for performing simulations of systems of particles.
The code is called 'lennyffs' since it was originally designed to do
forward flux sampling (FFS) simulations of systems with Lennard-Jones
(LJ) potentials.  However, the current capabilities of the code are
rather greater.  In particular, the code can be used with the
following potentials and algorithms:

potentials:
* Lennard-Jones
* Gaussian (also known as 'Gaussian Core')
* Inverse Power Law (IPL)

algorithms:
* NVT Metropolis Monte-Carlo
* NPT Metropolis Monte-Carlo
* NVE Molecular Dynamics

For rare event simulation, either FFS or Umbrella Sampling (US) can be
used,but 'direct simulation' (i.e. normal Monte-Carlo simulation) is
also possible.

REQUIREMENTS
------------

Python 2.x

Numpy

Boost C++ libraries

C++ compiler

Fortran compiler

INSTALLATION
------------

To install, simply run the script setup.sh

Below and in the rest of the document [base] refers to the base path
of the root code directory (the directory that this README.md file and
the setup.sh file is contained in).

Running setup.sh will do the following:

* Navigate to [base]/modules/fortran and run 'make' to compile fortran library

* Navigate to [base]/modules/cpp and run 'make' to compile C++ library

* Add the following paths to your PYTHONPATH environment variable:
  * [base]/modules/python
  * [base]/modules/fortran
  * [base]/modules/cpp

* Add the following directories to your PATH
  environment variable:
  * [base]/exec - this means that the executable files used
                such as [base]/exec/code.py can be run directly
					 from the command line with 'code.py' rather than 
					 having to type the full path.
  * [base]/util - same idea as above, but for accessing the utilities, 
                such as 'len_largestcluster' etc. directly from the 
					 command line.

CODE TREE
---------

* [base]/exec

* [base]/modules
  * /python
    * /ase
  * /cpp
    * /ops
  * /fortran
    * /clist
    * /gauss
    * /global
    * /len
    * /ops
    * /util

* [base]/other

* [base]/scripts

* [base]/util

RUNNING CODE
-------------

The following scripts are in the [base]/exec directory, marked executable (the first line is 
'#! /usr/bin/env python'):

code.py        - for direct simulations

codeffs.py     - for FFS simulations.  Note that codeffs.py submits jobs to a compute cluster.  It therefore must be run on a cluster computer with the Sun Grid Engine (SGE) installed (!), or a compatible engine that provides the 'qsub' command.

ffsdiagnostics.py - run after codeffs to collect FFS results (this can be run on any machine, no 'qsub' job submission needed).

usdiagnostics.py  - run after wumbrella.py to output opval (histogram) files.

NOTE: A number of sample input files are provided in [base]/examples
which illustrate how to run the code for both 'direct' and FFS
simulations.


code.py

This will look for an input file named 'in' in the working directory
it is run from (this should usually be a directory set up specifically
for the simulation, containing only this 'in' file).  The input file
contains details such as: whether the simulation should start afresh
or reload particles positions (and velocities if we are doing an MD
simulation) from a file, the potential to use (Lennard-Jones or
Gaussian), the number of particles, the box size, the temperature, the
simulation algorithm to use, the order parameter to use, how often
positions should be written to disk, etc.

codeffs.py

Again this will look for the input file name 'in' in the working
directory it is run from.  Note that this script CAN ONLY BE RUN FROM
A SYSTEM WHICH IS RUNNING THE SUN GRID ENGINE (typically the front-end
of a cluster computer which is responsible for handling the job
submission).  It uses the command 'qsub' to send multiple jobs to the
cluster.  All of the jobs for the entire FFS simulation (typically a
lot of jobs) are sent; codeffs.py uses the fact that some jobs can
'hold' until other jobs have completed, which means that e.g. the FFS
'shots' for interface 2 do not start before the interface 1 -> 2 has
finished.  See the file codeffs.py itself for further details.

umbrella.py

Again this will look for the input file name 'in' in the working
directory it is run from. Works in the same way as code.py except
the system can be subjected to a biasing potential to have some
control over which region of configuration phase can be explored.
While this can be run directly, it is usually more useful to run
wumbrell.py.

wumbrella.py

Will set up a series of umbrella sampling windows based on the
'in' file contents. It will also write a bash script to be
executed on the head node of a cluster running the Sun Grid Engine.
This script will run each of the windows simultaneously using
umbrella.py.

ffsdiagnostics.py

This should be run from a folder that codeffs.py has already been run
in, at the end of the simlation.  It will read the 'interface*.out'
files, which are written by codeffs.py at each interface, and collate
the statistics, writing out the number of successful shots at each
interface, and other diagnostics such as the transition rate.  Two
files are written, 'allinterfaces.out' and 'allinterfaceseff.out'.

usdiagnostics.py

This should be run from a folder in which a US simulation has been run
(typically the US simulation will have been run via wumbrella.py).  It
will read each of the opval*.out files (one for each window) written
by the US simulation and write fhist*.out files, which give the
(relative) free energies for the particular window.

INPUT FILE
----------

A number of sample input files are provided in [base]/examples which
illustrate how to run the code for both 'direct' and FFS simulations.
Each non-empty line of the input file that does not start with the '#'
character is read in.  The line should be of the form:

parameter value

Where parameter is the name of the input parameter, and value is its
value (which must be able to be converted to the correct type).  For
example, here are some possible examples of lines that appear in the
input file:

restartfile config1.xyz

potential len

Tstar 0.5

surface False

A complete list of the possible parameters, as well as their 'types'
(INT, FLOAT, STRING etc) can be found in
[base]/modules/python/paramdata.py .  Note that Python will raise
ValueError if the value of the parameter cannot be converted to the
type listed.

The combination of input parameters needed to run any given simulation
is slightly complicated and not well documented at the moment.  The
best thing to do is to look at the example input files in
[base]/examples .  After that, the details (source code) of how the
parameters in the input file are parsed are all contained in
[base]/modules/python/initsim.py .
