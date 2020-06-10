DESCRIPTION:
============

(N)on(E)quilibrium (A)nalysis of (R)ibo-seq is a program that finds codon-dependent elongation rates 
relative to the initiation rate from A-site density measured by Ribo-seq. 

MODEL:
======

The totally asymmetric simple exclusion process (TASEP) proposed by MacDonald CT, Gibbs JH, Pipkin AC. 
Kinetics of biopolymerization on nucleic acid templates. Biopolymers. 1968 Jan;6(1):1–25.

METHOD:
=======

The program uses nonlinear optimisation (least-squares) to find elongation-to-initiation rates that 
minimise the following objective function:

  S = [rho_2(exp) - rho_2(TASEP))]^2 + ... + [rho_L(exp) - rho_L(TASEP))]^2

where 

  rho_i(exp) = experimental A-site density at codon i
  rho_i(TASEP) = predicted A-site density at codon i

FILES:
======

/output       --> results
/ribo-seq     --> Ribo-seq data 
/source       --> Fortran source files
genelist.dat  --> list of gene names and average ribosome densities for normalisation
NEAR.dat      --> input parameters for NEAR
NEAR          --> executable program (for Linux)

File genelist.dat must have the following format:

  <genename> <average ribosome density>
  
where 
 
  <genename> = name of the gene (maximum 10 characters)
  <average ribosome density> = number of ribosomes / number of codons (not including START codon, including STOP codon)
   
File NEAR.dat has the following format

  <ll>			
  <mftol>		
  <psmtol>		
  <ftol>		
  <maxtime>
  <iter>
  <inputfile>

where

  <ll> = ribosome footprint lenght (default = 10)
  <mftol> = tolerance for the relative error between experimental density and theoretical density 
            predicted by the mean-field approximation (default = 0.05)
  <psmtol> = tolerance for the relative error between theoretical density predicted by the mean-field 
             approximation and by the power-series method (default = 0.1)
  <ftol> = stop iterating when the relative change in the objective function is at most ftol (default = 1.0d-8)
  <maxtime> = maximum time in seconds for optimisation of one gene (default = 3600)
  <iter> = number of iterations for Gillespie algorithm (default = 10000)
  <inputfile> = file with a list of genes to optimise (default = genelist.dat) 
  
RIBO-SEQ DATA
=============
 
Ribo-seq data is in folder /ribo-seq. The data must be in the following format:
 
  <codon> <number of A-site reads on that codon>

where 

  <codon> = codon (max 3 characters)
  <number of A-site reads on that codon> = number of A-site reads for that codon (must be integer)
  
There must be one file for each gene. Filename must be of the format "A-site-<genename>.dat. The same 
<genename> must be used in the genelist.dat file.
  

HOW TO COMPILE FROM SOURCE
==========================

In order to compile the program in Fortran, NLopt library is required (https://nlopt.readthedocs.io/en/latest/).
After installing NLopt, compile the Fortran file NEAR.f90 in the source directory following these steps:

  1. export LD_LIBRARY_PATH=<absolute path to nlopt>nlopt/lib64
  2. gfortran mt19937.f90
  3. gfortran NEAR.f90 -I<absolute path to nlopt>/nlopt/include -L<absolute path to nlopt>/nlopt/lib64 -lnlopt -lm -o NEAR mt19937.o
  4. mv NEAR ../.

Make sure that you copy the NEAR exectuable file to the main directory (step 4 above). NEAR file must be in the main directory.

HOW TO RUN THE PROGRAM
======================

In linux: ./NEAR from the main directory.