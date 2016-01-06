# Welcome to the PBR Tier System Code Repository
This repository currently includes Fortran 90/95 code for an age- and sex-structured operating model. The operating model can mimick the dynamics of either a single stock, or two stocks of the same species which overlap spatially. This allows for a range of management strategy evaluations to be performed. 

The code files are located in the `PBR Netbeans` folder. <a href="https://github.com/John-Brandon/PBR-Tier-System/tree/master/PBR%20Netbeans" target="_blank">The Readme.md file in that folder</a> contains descriptions for the code files. Additionally, R code is provided, as an example work-flow, for running trials as part of a management strategy evaluation for a single stock PBR Tier System.  

The Fortran code has been developed under Mac OS 10.9.5, using <a href="https://gcc.gnu.org/wiki/GFortran" target="_blank">the free, open source, GNU Fortran 95 (gfortran) compiler</a>. Instructions are provided below for compiling the Fortran code using the gfortran compiler. 

* For instructions on compiling the Fortran code under Mac OS X, [please see below](#compiling-the-fortran-code-under-mac-os)

* For instructions on compiling the Fortran code under Windows, [please see below](#compiling-the-fortran-code-under-windows)

Links are provided at the top of this page to view the files in this repository, as well as to download the files therein. 

### Funding and Contributors
This project is funded by the [U.S. Western Pacific Fisheries Management Council](http://www.wpcouncil.org/about-us/). Drs. John R. Brandon, [Andr&eacute; E. Punt](http://fish.washington.edu/people/punt/index.html), Paula Moreno and Randall R. Reeves are collaborators. They are working together as [the Independent Advisory Team on marine mammal stock assessment] (http://www.usm.edu/gcrl/scemfis/marine.mammal.asssessment.team.php), with the goal of developing research projects that serve to reduce uncertainty in marine mammal stock assessment. The team is funded through [Science Center for Marine Fisheries](http://scemfis.org/aboutus.html). The center is run under the [U.S. National Science Foundation's Industry and University Cooperative Research Program.](http://www.nsf.gov/eng/iip/iucrc/program.jsp) 


## Compiling the Fortran code under Mac OS X
(1). At least some versions of Mac OS X come with the gfortran (v4.2) compiler pre-installed. You can check if your machine has this installed by opening the terminal and typing
```shell
gfortran -v
``` 
(2). If you see a message like the one below, your machine has gfortran installed, and you can skip to step 5:
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
(4). Then re-type the `gfortran` command. If that doesn't work, you'll need to install gfortran. Links and instructions can be found at: http://hpc.sourceforge.net/

(5). Assuming you have gfortran installed, the next step is to make sure that you are in the working directory with the Fortran files downloaded from this repository. You can use the `cd` terminal command to change directories if you need.
```shell
cd Type-your-working-directory-here
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
GNU Fortran (GCC) 4.7.0 20111220 (experimental)
Copyright (C) 2011 Free Sofware Foundation, Inc.
```
(3). If you get an error message (e.g. the `gfortran` command can not be found), download MinGW from: http://sourceforge.net/projects/mingw

(4). Follow the installation steps in this video: https://www.youtube.com/watch?v=oVfAU1ziOjg

(5). Assuming you have gfortran installed, the next step is to make sure that you are in the working directory with the Fortran files downloaded from this repository. You can use the `cd` terminal command to change directories if you need.
```shell
cd Type-your-working-directory-here
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



