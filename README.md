PBR-Tier-System
===============

# General notes
This repository includes Fortran 90/95 code for a two-stock operating model. The operating model can also be a single stock by changing the `n_stocks` variable in the `input.par` file. 

The code is being developed under Mac OS 10.9.5, using [the *gfortran* 4.2 compiler (free, open source, GNU Fortran 95 compiler)](https://gcc.gnu.org/wiki/GFortran); it is being developed using the Netbeans 8.0 integrated development environment (IDE).   

The Fortran code is currently located in the `PBR Netbeans` folder of this repository.

# Status of code development
The operating model is nearly complete. I'm testing the code to allocate spatial bycatch (the population dynamics and Tier 2 (Wade, 1998) survey estiamtes of abundance have been tested). There is also R code being developed for reading output files and creating plots, but I haven't added that code to this repository (yet)...



