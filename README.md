# Welcome to the PBR Tier System Code Repository
This repository currently includes Fortran 90/95 code for an age-structured multi-stock population dynamics model. Initial code development is occurring, under Mac OS 10.9.5, using <a href="https://gcc.gnu.org/wiki/GFortran" target="_blank">the free, open source, GNU Fortran 95 compiler</a>. In the near future, a related open-source repository will be linked with R code to visualize and analyze the output files from the Fortran code. 

Links are provided at the top of this page to view the files in this repository, as well as to download the files therein (raw code or otherwise). Motivation and background for this project is provided below, but first a quick FAQ.

## FAQ
__Q1:__ 

What does "open source" meen?

__A1:__ 

You have the power to fix the typo in __Q1__. 

All of the files in this repository (Fortran, R, markdown, or otherwise) are available to the public, including the markdown text file underlying this webpage. Git provides us a formal process of version control for collaborating on these files. 

If you're new to Git and version control, and interesting in learning how to collaborate on GitHub, <a href="https://www.youtube.com/watch?v=SCZF6I-Rc4I" target="_blank"> here's a nice video (6min 33s) with an example tututorial. </a> As a first exercise, you could access  <a href=https://github.com/John-Brandon/PBR-Tier-System/tree/master target="_blank"> the markdown text for this webpage in the README.md file </a>, following the instructions in the video. Once you've edited the typo above, submit your first pull request. There's no better way to learn than by doing. 

Regardless of whether you know how to use Git, if you see something that is not factually correct; important information is missing; find a bug in the code; have a better way of doing something -- let us know (<a href="https://guides.github.com/features/issues/" target="_blank">here's one avenue through GitHub</a>). You can also contact the maintainer of this repository, listed at the end of this page. 

Legally speaking, for this repository, open source means that the files here are licensed under the <a href="https://github.com/John-Brandon/PBR-Tier-System/blob/master/LICENSE.md" target="_blank">MIT License</a>. You can use them however you want, even for profit, but you can't hold us liable for any damages. 

__Q2__: 

Why is open source important?

__A2__: 

One of the goals of this project is to provide an opportunity to engage and grow the scientific community for simulation research on "Potential Biological Removal" (*aka* PBR -- More on PBR below). We've decided to use GitHub as an open-source platform to ensure that this research is transparent and reproducible. 

We hope this project excites ideas and methodological approaches for reducing scientific uncertainty in the context of meeting the management objectives of the U.S. Marine Mammal Protection Act, or other related conservation objectives. We also offer this as a platform to test existing strategies to calculate PBR, which haven't been formally tested yet. 

PBR is also a useful concept for scientists outside the U.S., who are looking to provide management advice on sustainable levels of human caused mortality. An open-source repository engenders international collaboration and improvements in scientific advice. 

Finally, you just never know who might be reading this and want to collaborate on a scientifically interesting and challenging problem. For all we know, the reader might be a curious computer scientist working for a tech company, with a desire learn more about computational applications in natural resource science. Making this project open source allows anyone with an internet connection to contribute to the code and development of algorithms.

With all of that said, let's get into the fun, starting with some background...

## The U.S. Marine Mammal Protection Act and Potential Biological Removal
The [U.S. Marine Mammal Protection Act (MMPA)](http://www.nmfs.noaa.gov/pr/laws/mmpa/text.htm) contains rules to limit the number of human caused mortalities to individual stocks of marine mammals (stocks are the unit to be managed). In particular, the legislation prescribes a "Potential Biological Removal" limit to human caused mortality. Further, the MMPA mandates that this limit should be calculated as the product of three values: (1) An estimate of abundance that "provides reasonable assurance that the stock size is equal to or greater than the estimate"; (2) One half of the maximum intrinsic population growth rate; and (3) A recovery factor between 0.1 and 1.0. 

PBR can therefore be written as the equation: 
PBR = N<sub>MIN</sub> 0.50 R<sub>MAX</sub> F<sub>R</sub> 

The management objective of PBR, as outlined in the MMPA and quantified by the U.S. National Marine Fisheries Science Service (NMFS), is to allow stocks of marine mammals to be maintained at an abundance that is 50% or more of the carrying capacity of their environment. 

## Previous simulation research on catch control rules 
Simulation testing of catch control rules was introduced to fisheries science by the Scientific Committee of the International Whaling Commission (the SC). The SC's research and development started during the 1980s and is ongoing ([see Punt and Donavon, 2007, for a review](http://icesjms.oxfordjournals.org/content/64/4/603.full.pdf+html)). 

Generally speaking, the procedures that have been developed take into account scientific uncertainties as represented by a range of 'what if' scenarios. For example, "What if abundance estimates are biased?" Would the catch control rule still meet the management objectives? Each scenario is known as a "trial". Simulation modeling is used to evaluate the performance of alternative management strategies (e.g. catch control rules), across a set of trials, given quantified management objectives. 

This approach for providing scientific advice to fisheries management has also become increasingly popular in recent years, and is now also used for evaluating management systems of commercial fisheries for fin- and shell-fish. The procedures originally developed by the SC are now known in fisheries science as "Management Strategy Evaluation" (MSE). 

## Management strategy evaluation 
MSE involves a few key steps, including: (i) development of an "operating model", which represents the system being managed (the 'true' state of nature) and how observable data are generated; (ii) identification of a candidate management strategy, which also includes a data collection strategy;  (iii) developing a set of trials that span the range of key uncertainties about the true state of nature, including trials where data are biased and/or imprecise; and (iv) evaluation, using simulation, of the candidate management strategies. Candidate management strategies are evaluated based on their performance relative to the quantified management objectives (e.g., maintaining abundance above some level with 95% certainty). In the context of this project, the management objectives of the U.S. Marine Mammal Protection Act are the benchmark. 

## Previous simulation research on PBR
Initial work by NMFS scientists on simulation testing of PBR was born from the MSE procedures developed by the SC. NMFS scientists are also active members of the SC. 

[Taylor (1993)](https://swfsc.noaa.gov/publications/TM/SWFSC/NOAA-TM-NMFS-SWFSC-188.PDF) used simulation to explore the impact of assigning to N<sub>MIN</sub> alternative percentiles of an abundance estimate, given survey sampling error (e.g. assigning the mean vs lower 95th confidence interval of the abundance estimate to N<sub>MIN</sub>). Taylor (1993) demonstrated that a strategy of assigning the point estimate of abundance to N<sub>MIN</sub> resulted in poor performance when biase existed in the abundance data. A strategy of assigning the lower 95th confidence interval of the abundance estimate to N<sub>MIN</sub> was relatively robust to the biases that were explored. 

Wade (1998) expanded on the simulation studies by Taylor (1993). As a first step, Wade (1998) iterated over a set of trials where estimates of abundance were unbiased, but in some cases imprecise. For these trials, Wade (1998) essentially solved for the percentile of the abundance estimate that resulted in meeting the management objective of the MMPA in 95% of the simulations across scenarios. That value was the lower 20th percentile of the abundance estimate. Wade (1998) then examined a variety of trials where different biases existed in the data. For these trials, N<sub>MIN</sub> was set to the lower 20th percentile of the abundance estimate, and the process was repeated to find the solution for the value of F<sub>R</sub> that resulted in 95% of the simulations meeting the management objectives if data are biased. The results of Wade's simulation testing form the basis for the values NMFS currently uses when calculating PBR for stocks of marine mammals. 

Below, we describe how the PBR Tier System study aims to expand on previous research. We also provide some rationale for doing so. 

TODO: Add some general rationale for a tier system, before getting into modeling. Could link the poker example here: https://github.com/John-Brandon/PBR-Tier-System/blob/gh-pages/poker.md

## Operating model synopsis
The underlying population dynamics were represented using an age-, sex-, and space-structured multiple-stock model (which also allows single-stock scenarios to be evaluated). The age-structured component of the model is an extension of single-stock assessment models employed by the NMFS scientists and the SC of the IWC (e.g., Punt 1999, Wade 2002). The single-stock model is expanded here to allow for two spatially overlapping stocks of the same species (Fig 1). This approach is similar to multi-stock models developed by the SC of the IWC for MSEs under more complicated scenarios, involving overlapping stocks, and uncertainties in stock identity with respect to catch limits (e.g. IWC 2012).

![Spatial Structure of the Model](https://dl.dropboxusercontent.com/u/31861309/SpatialFig.png)

*Figure 1.* The ranges and overlap of the two stocks are shown over the four areas. The surveyed area comprises areas 1, 2, and 3. Estimates of abundance include all animals (both stocks) in the survey area. The ranges of the two stocks overlap in area 2. 

The modeling approach differs from the age-aggregated model used in previous PBR research (e.g., Taylor 1993; Wade 1998). For simple scenarios, the two approaches should result in equivalent dynamics. However, the more complex population dynamics model developed here allows a wider spectrum of trials to be evaluated (e.g., spatial overlap of stocks of the same species, and non-uniform age- and sex-specific human caused mortality).

The Tier System aspect of this research involves different levels of data availability. Previous PBR research has largely  assumed that PBR would be calculated based on a single absolute estimate of abundance. 

(although see Wade and DeMaster 1998, for an example of a PBR simulation study where multiple estimates were averaged: TODO Add this reference and wrap up the remaining last bits of text)

TODO?: Add image from slides showing hierarchy of tiers. 
(in this case, features of the proposed tier system might include how historical abundance estimates are weighted, how trends are estimated, and whether abundance data older than 8 years are used)

TODO: Add preliminary image showing population trajectories. 

### Funding and Contributors
This project is funded by the [U.S. Western Pacific Fisheries Management Council](http://www.wpcouncil.org/about-us/). Drs. John R. Brandon, [Andr&eacute; E. Punt](http://fish.washington.edu/people/punt/index.html), Paula Moreno and Randall R. Reeves are collaborators. They are working together as [the Independent Advisory Team on marine mammal stock assessment] (http://www.usm.edu/gcrl/scemfis/marine.mammal.asssessment.team.php), with the goal of developing research projects that serve to reduce uncertainty in marine mammal stock assessment. The team is funded through [Science Center for Marine Fisheries](http://scemfis.org/aboutus.html). The center is run under the [U.S. National Science Foundation's Industry and University Cooperative Research Program.](http://www.nsf.gov/eng/iip/iucrc/program.jsp) 

### References 
IWC. 2012. International Whaling Commission. 2012. Report of the Scientific Committee. Annex E. Report of the Standing Working Group on an Aboriginal Subsistence Whaling Management Procedure. J. Cetacean Res. Manage. (Suppl.) 13:130-53. 

Punt, A. E. (1999). A full description of the standard Baleen II model and some variants thereof. J. Cet. Res. Manage. 1:267-276

[Taylor, B.L. 1993. “Best” abundance estimates and best management: Why they are not the same. U.S. Department of Commerce, NOAA Technical Memorandum NMFS-SWFSC-188. 20 pp.](https://swfsc.noaa.gov/publications/TM/SWFSC/NOAA-TM-NMFS-SWFSC-188.PDF)

Wade, P.R. 1998. Calculating limits to the allowable human-caused mortality of cetaceans and pinnipeds. Mar Mamm Sci. 14: 1-37.

Wade, P.R. 2002. A Bayesian stock assessment of the eastern Pacific gray whale using abundance and harvest data from 1967-1996. J. Cet. Res. Manage. 4: 85-98. 

---------------------------

*last updated: 10 Jun 2015*

*contact: jbrandon@gmail.com*
