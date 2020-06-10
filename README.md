# NEAR
Non-Equilibrium Analysis of Ribosome profiling data

## Description

Non-Equilibrium Analysis of Ribo-seq (NEAR) is a program that infers codon-dependent elongation rates relative to the initiation rate from A-site density measured by ribosome profiling. The program uses nonlinear optimisation to match the experimental ribosome density to the theoretical density predicted by a mathematical model for ribosome dynamics called the totally asymmetric simple exclusion process (TASEP) \[1\]. The progam is described in detail in \[2\].

## Installation

The program is written in Fortran. Before compiling the program, [NLopt library](https://nlopt.readthedocs.io/en/latest/) must be installed. NLopt library is available for Unix-like systems and Windows. After installing NLopt, the program can be compiled using `gfortran`. On a Linux system, go to the directory `source` and set the path to the NLopt library:

> export LD_LIBRARY_PATH=\<*absolute path to NLopt*\>/lib64
  
where \<*absolute path to NLopt*\> is the absolute path to the NLopt directory. Then compile the object file for the Mersenne Twister pseudorandom number generator:

> gfortran -c mt19937.f90

Finally, compile the program named NEAR linking to the NLopt library:

> gfortran NEAR.f90 -I\<*absolute path to nlopt*\>/include -L\<*absolute path to nlopt*\>/lib64 -lnlopt -lm -o NEAR mt19937.o

After compiling, move the executable file NEAR to the main directory.
  
## Input data

### NEAR.dat

The main parameters of the program are stored in the file *NEAR.dat*. The parameters are:

- `l`: ribosome footprint lenght (*default* = 10)
- `mftol`: tolerance for the relative error between experimental density and theoretical density predicted by the mean-field approximation (*default* = 0.05)
- `psmtol`: tolerance for the relative error between theoretical density predicted by the mean-field approximation and by the power-series method (*default* = 0.1)
- `ftol`: stop iterating when the relative change in the objective function is at most ftol (*default* = 1.0d-8)
- `maxtime`: maximum time in seconds for optimisation of one gene (*default* = 3600)
- `iter`: number of iterations for Gillespie algorithm (*default* = 10000)
- `inputfile`: file with a list of genes to optimise (*default* = genelist.dat)

### genelist.dat

The file *genelist.dat* contains names of genes to optimise. The format of the file is:

> \<*genename*\> \<*average ribosome density*\>    
> \<*genename*\> \<*average ribosome density*\>  
> ...    

where \<*genename*\> is a name for the gene (maximum 10 characters) and \<*average ribosome density*\> is equal to the number of ribosomes per transcript length, where the transcript length is equal to the number of codons including STOP but exclusing START.

### ribo-seq directory

*ribo-seq* directory should contain files with A-site footprint reads. Each file should be named *A-site_\<genename\>.dat*, where \<*genename*\> should be matched to the \<*genename*\> in the file *genelist.dat*. The format of the file *A-site_\<genename\>.dat* is:

> \<*codon*\> \<*number of A-site footprint reads mapped to that codon*\>      
> \<*codon*\> \<*number of A-site footprint reads mapped to that codon*\>       
> ...        

## Output data

All output data are written in the directory *output*. For each gene in the *genelist.dat* file, two output files are generated. One file are the results of the optimisation procedure, called \<*genename\>-results.dat*. The other file is a log file, where details of the optimisation procedure are stored, for example if the optimisation was successful. Results are stored in the following format:

> \<*codon*\> \<*included*\> \<*rhoexp*\> \<*rhomc0*\> \<*rhomc*\> \<*rhopsm*\> \<*current0*\> \<*current*\> \<*f1*\> \<*f2*\> \<*f3*\> 

where

- `codon`: codon identity
- `included`: 0 if the codon was not further optimised (meaning that the initial guess was good enough), 1 if it was optimised
- `rhoexp`: experimental density normalised using \<*average ribosome density*\> from *genelist.dat* file
- `rhomc0`: theoretical density obtained by Monte Carlo simulations using the initial rates 
- `rhomc`: theoretical density obtained by Monte Carlo simulation using the optimised rates
- `rhopsm`: theoretical density obtained by power series method using the optimised rates
- `omega0`: initial elongation-to-initiation rate 
- `omega`: optimised elongation-to-initiation rate at codon `codon`
- `current0`: predicted translation rate relative to the initiation rate obtained using the initial rates 
- `current`: predicted translation rate relative to the initiation rate obtained using the optimised rates
- `f1`: first-order coefficient for the power series expansion of the ribosome density
- `f2`: second-order coefficient for the power series expansion of the ribosome density
- `f3`: third-order coefficient for the power series expansion of the ribosome density

The main results are `omega`, `rhomc` and `current`. Other results are for the quality check of the optimisiation procedure. See \[2\] for more details.


## References

\[1\] MacDonald CT, Gibbs JH, Pipkin AC. Kinetics of biopolymerization on nucleic acid templates. [Biopolymers. 1968 Jan;6(1):1–25.](https://doi.org/10.1002/bip.1968.360060102)
\[2\] Szavits-Nossan J, Ciandrini L. Inferring efficiency of translation initiation and elongation from ribosome profiling. [bioRxiv](https://doi.org/10.1101/719302)


