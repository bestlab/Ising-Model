# Ising-like Model
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

Example conf file and villin contact data found in example folder.
Change to example and run ../bin/run_ising example/villin.conf
The output should be similar to the villin.bm_param file found in Zenodo.

To run a continuation simulation use a conf file similar to the villin_cont.conf file found in /example.
The main differences between a conf file for running the initial simulation and a conf file for running a 
continuation simulation:
	
 	1. freq_dir: uses frequency files found in scratch directory of initial simulation, specifies
 	the frequencies of single and pairwise contacts found in the msa_name file so it does not have
  	to be recalculated
   	2. input_name: uses previous parameter files found in scratch directory of initial simulation,
	specifies the parameters used in the last iteration of the previous simulation
	3. mc_steps: this should be updated to be the same as the maximum number of mc steps used
	in the previous simulation. This information can be found in the log file under lines detailing
	"MC steps per iteration: "
	4. num_restarts: update to indicate number of continuations

Background:

The run_ising program uses the maximum likelihood principle to infer a probability
distribution for the energetics of contact formation. This probability
distribution is in the form of an Ising Model and is characterized by an energy
function consisting of single contact energy and two-site contact-contact pair
energy terms. The code uses Boltzmann machine learning (following Lapedes,
Giraud, & Jarzynski arXiv:1207.2484) to learn these parameters given a set of
trajectory frames where all possible contacts have been delineated in the code.
To run the code, build the code and run the executable in the bin folder. 
The script requires as input a conf file; one example of which is located in 
the "example" folder. The conf file specifies input and output directories 
and the location of the MD trajectory for the input data.

All simulation code is written in C++ and has been successfully compiled by the
author with C++17 and C++14 compilers. The code may not compile on every
computer setup, may require modifications to source code or installation of
additional libraries/compilers to compile, and may not be robust to all input
parameters. The ising code has been tested for its ability to produce energy
models which, when sampled, can reproduce the single and pairwise contact frequencies
of an input trajectory to within a user-defined accuracy. However, for very high
accuracy (low error tolerance), the code may take infeasibly long to run. This
may be especially true for proteins with longer sequences. Furthermore, given
the quadratic scaling of the number of contacts with respect to the number of
residues in a protein we estimate that the code only functions for proteins
with ~180 residues or less. Note that compilation (which can be achieved with 
the "make" command in the main directory) and use require that the 
Armadillo linear algebra library for C++ be installed.

