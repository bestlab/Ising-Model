# Boltzmann-Machine-DCA
Implementation of ising-like model for protein contact analysis using a Boltzmann machine algorithm.

There is one executable produced in this repository: run_ising. See below for background.

Installation Prerequisites:
	1. C++ compiler and standard libraries
	2. GNU Scientific Library (GSL) https://www.gnu.org/software/gsl/ 
	3. Armadillo C++ library for lineary algebra & scientific computing https://arma.sourceforge.net/

Installation Instructions:
	1. Edit makefile to set local paths, in particular:
		- ARMA_DIR set to installation of Armadillo library
		- CXXFLAGS and LDFLAGS set to include any nonstandard paths for libraries
	2. Run make to compile

Test Example Instructions:

Change to example and run ../bin/run_ising example/ntl9_stable_eg.conf

Analysis scripts:

coop_distribution_tails_log.py
	-- computes distribution of cooperativities from output of Ising Model fit

coop_network.py
	-- computes cooperative contact networks from output of Ising Model fit

Background:

The run_ising program uses the maximum likelihood principle to infer a probability
distribution for the energetics of contact formation. This probability
distribution is in the form of an Ising Model and is characterized by an energy
function consisting of single contact energy and two-site contact-contact pair
energy terms. The code uses Boltzmann machine learning (following Lapedes,
Giraud, & Jarzynski arXiv:1207.2484) to learn these parameters given a set of
trajectory frames where all possible contacts have been delineated in the code.
To run the code, use the "run_ising.sh" script located in the
"Ising-Model/scripts" folder. The script requires as input a conf file, several
examples of which are located in the "swDCA" folder. The conf file specifies
input and output directories and files as well as learning parameters. An MSA
file must be specified and provided. Several examples are located in the
"swDCA" folder. Currently, the code supports MSAs in FASTA (.fas) format, as
well as MSAs of HP model sequences (.hp). 

All simulation code is written in C++ and has been successfully compiled by the
author with C++17 and C++14 compilers. The code may not compile on every
computer setup, may require modifications to source code or installation of
additional libraries/compilers to compile, and may not be robust to all input
parameters. The swDCA code has been tested for its ability to produce energy
models which, when sampled, can reproduce the one- and two-site frequencies of
the input MSA to within a user-defined accuracy. However, for very high
accuracy (low error tolerance), the code may take infeasibly long to run. This
may be especially true for longer protein sequences. Note that compilation
(which can be achieved with the "make" command in the "swDCA" directories) and
use require that the Armadillo linear algebra library for C++ be installed.

Analysis scripts (located in the "scripts" subdirectories of "swDCA") are
written in Python 3 and bash. 

