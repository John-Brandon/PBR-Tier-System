# Welcome to the PBR Tier System Code Repository
This repository currently includes Fortran 90/95 code for an age- and sex-structured operating model. Additionally, R code is provided which controls the batching of simulation trials, compiles performance statistics, and plots output.

The operating model can mimick the dynamics of either a single stock, or two stocks of the same species which overlap spatially. A full mathematical description of the operating model <a href="https://www.dropbox.com/sh/qga2x5sq2h41vfp/AADlfFXYeO9MjfjrRor4-Z1Ca?dl=0" target="_blank">is available here</a>

This code allows for a tier system approach to management strategy evaluation. The PBR formula as evaluated by Wade (1998) makes use of only the most recent estimate of abundance, irrespective of its precision and the number of estimates available. A tier system would assign species or stocks to different tiers based on the availability and quality of data. Ultimately, a tier system approach would serve to optimize available information for each stock.

The Fortran and R code files are located in the `PBR Netbeans` folder. <a href="https://github.com/John-Brandon/PBR-Tier-System/tree/master/PBR%20Netbeans" target="_blank">The Readme.md file in that folder</a> contains an overview description of each code file.  

The Fortran code has been developed under Mac OS 10.9.5, using <a href="https://gcc.gnu.org/wiki/GFortran" target="_blank">the free, open source, GNU Fortran 95 (gfortran) compiler</a>. Instructions are provided below for compiling the Fortran code using the gfortran compiler. Compiled executables for Mac OS X and Windows <a href="https://www.dropbox.com/sh/qga2x5sq2h41vfp/AADlfFXYeO9MjfjrRor4-Z1Ca?dl=0" target="_blank">are available here.</a> 

### Funding and Contributors
This project was funded by the <a href="http://www.wpcouncil.org/about-us/" target="_blank">U.S. Western Pacific Fisheries Management Council</a>. Drs. <a href="https://www.linkedin.com/in/john-brandon-b5690a26" target="_blank">John R. Brandon</a>, <a href="https://fish.uw.edu/faculty/andre-punt/" target="_blank">Andr&eacute; E. Punt, U.W.</a>, Paula Moreno (U.S.M), and Randall R. Reeves (Okapi Wildlife Assoc.) are collaborators. They are working together as <a href="http://scemfis.org/aboutus.html" target="_blank">Science Center for Marine Fisheries</a>, with the goal of developing research projects that serve to reduce uncertainty in marine mammal stock assessment. The team is funded through the <a href="http://www.nsf.gov/eng/iip/iucrc/program.jsp" target="_blank">text</a>. The center is run under the <a href="" target="_blank">U.S. National Science Foundation's Industry and University Cooperative Research Program.</a> 

## Compiling the Fortran code under Mac OS X
(1). Recent versions of Mac OS X may have the gfortran compiler pre-installed. You can check if your machine has this installed by opening the terminal and typing
```shell
gfortran -v
``` 
(2). If you see a message like the one below, your machine has gfortran installed, and you can skip to step 5 (noting step 3):
```shell
Using built-in specs.
Target: i686-apple-darwin8
Configured with: /Builds/unix/gcc/gcc-4.2/configure --prefix=/usr/local --mandir=/share/man --program-transform-name=/^[cg][^.-]*$/s/$/-4.2/ --build=i686-apple-darwin8 --host=i686-apple-darwin8 --target=i686-apple-darwin8 --enable-languages=fortran
Thread model: posix
gcc version 4.2.3
```
(3). If you get an error message (e.g. the `gfortran` command can not be found), try installing Apple's XCode command-line tools by typing:
```shell
xcode-select --install
```
(4). Then re-type the `gfortran -v` command. If that doesn't work, you'll need to install gfortran. Links and instructions can be found at: http://hpc.sourceforge.net/

(5). Assuming you have gfortran installed, the next step is to make sure that you are in the working directory with the Fortran files downloaded from this repository (e.g. ~/PBR-Tier-System-master/PBR Netbeans). You can use the `cd` terminal command to change directories.
```shell
cd Your-working-directory-here
# e.g.
# cd ~/PBR-Tier-System-master/PBR\ Netbeans
```
(6). Next, use this command to compile the Fortran code:
```shell
gfortran A_Random_module.f90 BRENT.f90 Declare_variables_module.f90 Eigen_module.f90 Generate_random_numbers_module.f90 Initialize_pop_module.f90 PBR_Errorcheck_module.f90 PBR_FileIO_Module.f90 PBR_calcs_module.f90 main.f90 -o main -fbounds-check -framework accelerate 
```
(7). Once everything compiles, you can run the program by typing:
```shell
./main
```

## Compiling the Fortran code under Windows
Windows doesn't come with gfortran pre-installed. Nevertheless, it's worth double checking that it hasn't been installed at some point in the past (perhaps bundled-in with another installation).  

(1). Open the command prompt and type: 
```shell
gfortran --version
```
(2). If you see a message like the one below, your machine has gfortran installed, and you can skip to step 5:
```shell
GNU Fortran (GCC) 4.7 20111220 (experimental)
Copyright (C) Free Sofware Foundation, Inc.
```
(3). If you get an error message (e.g. the `gfortran` command can not be found), download MinGW from: http://sourceforge.net/projects/mingw

(4). Follow the installation steps in this video: https://www.youtube.com/watch?v=oVfAU1ziOjg

(5). Assuming you have gfortran installed, the next step is to make sure that you are in the working directory with the Fortran files downloaded from this repository (e.g. C:\PBR-Tier-System-master\PBR Netbeans). You can use the `cd` terminal command to change directories if you need.
```shell
cd Type-your-working-directory-here
# e.g.
# cd C:\PBR-Tier-System-master\PBR Netbeans
```

(6). Next, use this command to compile the Fortran code:
```shell
gfortran A_Random_module.f90 BRENT.f90 Declare_variables_module.f90 Eigen_module.f90 Generate_random_numbers_module.f90 Initialize_pop_module.f90 PBR_Errorcheck_module.f90 PBR_FileIO_Module.f90 PBR_calcs_module.f90 main.f90 liblapack.a libblas.a -o main 
```

(7). Once everything compiles, you can run the program by typing:
```shell
main
```

[Top](#welcome-to-the-pbr-tier-system-code-repository)

*last updated: 2015 Jan 05*

*contact: jbrandon@gmail.com*



