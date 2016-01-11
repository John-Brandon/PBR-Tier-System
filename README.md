# Welcome to the PBR Tier System Code Repository
This repository currently includes Fortran 90/95 code for an age- and sex-structured operating model. Additionally, R code is provided which serves as a script to: automate the generation of input files for trials (<a href="http://www.cheatography.com/davechild/cheat-sheets/regular-expressions/" target="_blank">e.g. using regex</a>); batch multiple runs of simulation trials through shell commands; compile performance statistics across the batch, and; plot output. Instructions for installing R  <a href="https://cran.r-project.org/" target="_blank">are available here.</a> 

The operating model can mimick the dynamics of either a single stock, or two stocks of the same species which overlap spatially. A full mathematical description of the operating model <a href="https://www.dropbox.com/sh/qga2x5sq2h41vfp/AADlfFXYeO9MjfjrRor4-Z1Ca?dl=0" target="_blank">is available here</a>

This code allows for a tier system approach to management strategy evaluation. The PBR formula as evaluated by Wade (1998) makes use of only the most recent estimate of abundance, irrespective of its precision and the number of estimates available. A tier system would assign species or stocks to different tiers based on the availability and quality of data. Ultimately, a tier system approach would serve to optimize available information for each stock.

The Fortran and R code files are located in the `PBR Netbeans` folder. <a href="https://github.com/John-Brandon/PBR-Tier-System/tree/master/PBR%20Netbeans" target="_blank">The README.md file in that folder</a> contains an overview description of each code file and notes on running the R scripts.  

The Fortran code has been developed under Mac OS 10.9.5, using <a href="https://gcc.gnu.org/wiki/GFortran" target="_blank">the free, open source, GNU Fortran 95 (gfortran) compiler</a>. Instructions are provided below for compiling the Fortran code using the gfortran compiler. Compiled executables for Mac OS X and Windows <a href="https://www.dropbox.com/sh/qga2x5sq2h41vfp/AADlfFXYeO9MjfjrRor4-Z1Ca?dl=0" target="_blank">are available here.</a> 

### Funding and Contributors
Funding for this project was provided by the <a href="http://www.wpcouncil.org/about-us/" target="_blank">U.S. Western Pacific Fisheries Management Council</a>. Drs. <a href="https://www.linkedin.com/in/john-brandon-b5690a26" target="_blank">John R. Brandon</a>, <a href="https://fish.uw.edu/faculty/andre-punt/" target="_blank">Andr&eacute; E. Punt (U.W.)</a>, Paula Moreno (U.S.M), and Randall R. Reeves (Okapi Wildlife Assoc.) are collaborators, working together as part of the<a href="http://gcrl.usm.edu/scemfis/marine.mammal.asssessment.team.php" target="_blank"> Independent Advisory Team for Marine Mammal Assessment (IAT)</a>. The IAT was formed under the <a href="http://scemfis.org/aboutus.html" target="_blank">Science Center for Marine Fisheries</a>, with the goal of developing research projects that serve to reduce uncertainty in marine mammal stock assessment. The Science Center is managed under the <a href="http://www.nsf.gov/eng/iip/iucrc/program.jsp" target="_blank">U.S. National Science Foundation's Industry and University Cooperative Research Program.</a> 

## Compiling the Fortran code under Mac OS X
(1). Recent versions of Mac OS X may have the gfortran compiler pre-installed. You can check if your machine has this installed by opening the terminal and typing
```shell
gfortran --version
``` 
(2). If you see a message like the one below, your machine has gfortran installed, and you can skip to step 5 (noting step 3):
```shell
GNU Fortran (GCC) 4.2.3
Copyright (C) 2007 Free Software Foundation, Inc.

GNU Fortran comes with NO WARRANTY, to the extent permitted by law.
You may redistribute copies of GNU Fortran
under the terms of the GNU General Public License.
For more information about these matters, see the file named COPYING
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
## Compiling the Fortran code under Ubuntu Linux 
These instructions were tested on Ubuntu v14.04 (thanks to Michael Mathews for his trouble-shooting help with Linux). 

(1). Open the terminal and type: 
```shell
gfortran --version
```
(2). If you see a message like the one below, your machine has gfortran installed, and you can skip to step 5:
```shell
GNU Fortran (Ubuntu 4.8.4-2ubuntu1~14.04) 4.8.4
Copyright (C) 2013 Free Software Foundation, Inc.
```
(3). If you get an error message (e.g. the `gfortran` command can not be found) type:
```shell
sudo apt-get install gfortran
```
(4). Assuming you now have gfortran installed, the next step is to check whether you have the LAPACK and BLAS libraries installed. If you have R installed, LAPACK and BLAS should be installed already in the /usr/lib/ directory. You can check by typing:
```shell
locate liblapack.a
locate libblas.a
```
If you see those files in the usr/lib/ directory, you can move on to the next step. Otherwise, you'll need to install LAPACK and BLAS. These libraries are available from: http://www.netlib.org/lapack

(5). Assuming you have gfortran, LAPACK and BLAS installed, make sure you are in the working directory with the Fortran files downloaded from this repository (e.g. ~/PBR-Tier-System-master/PBR Netbeans). You can use the `cd` terminal command to change directories if you need.
```shell
cd Type-your-working-directory-here
# e.g.
# cd ~/PBR-Tier-System-master/PBR Netbeans
```
(6). Next, use this command to compile the Fortran code:
```shell
gfortran A_Random_module.f90 BRENT.f90 Declare_variables_module.f90 Eigen_module.f90 Generate_random_numbers_module.f90 Initialize_pop_module.f90 PBR_Errorcheck_module.f90 PBR_FileIO_Module.f90 PBR_calcs_module.f90 main.f90 -o main -L/usr/lib/ -llapack -lblas 
```
(7). Once everything compiles, you can run the program by typing:
```shell
./main
```

[Top](#welcome-to-the-pbr-tier-system-code-repository)

*last updated: 2015 Jan 05*

*contact: jbrandon@gmail.com*



