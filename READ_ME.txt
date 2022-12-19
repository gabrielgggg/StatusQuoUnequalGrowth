Code Appendix for 


  ``Bargaining over Taxes and Entitlements in the Era of Unequal Growth''
                Azzimonti, Karpuska, Mihalache (IER, 2022)


This packages includes
~~~~~~~~~~~~~~~~~~~~~~
[1] unequal.f90 : the Fortran source code for the main program
[2] Makefile : a sample compilation script/Makefile
[3] results/ : MATLAB scripts used to generate plots and further analysis
[4] twoPeriod/ : MATLAB code for the two-period version of our model


Compilation and use of the infinite horizon model code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main program for the infinite horizon version of the model is written in
Modern Fortran (free form) and parallelized using OpenMP. It was tested with:
 * Intel OneAPI Classic Compiler (ifort) version 2021.06 
    on Windows (with Visual Studio 2022) and on Linux (Ubuntu 22.04)
 * GNU gfortran 11 on Linux (Ubuntu 22.04)
 * Cray Compiler 10 for ARM on CentOS v8
 * Nvidia HPC SDK 22.7 on Windows and on CentOS v8
 
The MATLAB scripts for plotting were tested on R2021a and newer.

1. Compile the unequal.f90 Fortran source file, with OpenMP parallelization

2. Run the resulting executable
     The output of the code includes simulation statistics for the model.
     
3. The executable requires that a "results/" subfolder/subdirectory is
     present in the folder where the program is executed
     
4. Once completed, the code will write to disk, in results/, a .tab file
     containing tab-separated parameter values in a fixed order, and
     several .bin files containing the binary data for various values and
     policy functions.
     
5. In results/, the MATLAB script toMatlab.m loads these files and saves
     a MATLAB dataset containing all results, called results.mat
     (This script requires the loadBinary.m function.)
     
6. The transition.m script replicates the time paths based on turn over
     in US data, based on Senate, House, or President's Party as well
     as the case of ex-post persistent incumbency.
     
7. The samplePlots.m script provide examples of policy function plots

Alternative rules:
  On line 47 the code is configured to compute the equilibrium under one
  of the status quo rules defined on lines 41-46. Our baseline case is
  RULE_BOTH, where the status quo specifies both a tax level "tau" and an
  entitlement level "e". 

All other model parameters are set on lines 18-27. Numerical parameters---
governing convergence criteria, taste shocks, etc---are set on lines 30-34.


The two period model
~~~~~~~~~~~~~~~~~~~~

The MATLAB code for this model were tested on MATLAB R2011a and newer.

The main script for the two period model is ...TODO...





Contact
~~~~~~~

For question related to running this code, please contact the authors at:
  Marina Azzimonti: marina.azzimonti@gmail.com
  Laura Karpuska: laukarpuska@gmail.com
  Gabriel Mihalache: mihalache@gmail.com



