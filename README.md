# Welcome to the PBR Tier System Code Repository
<a href="https://github.com/John-Brandon/PBR-Tier-System/tree/master/PBR%20Netbeans" target="_blank">This repository currently includes Fortran 90/95 code for an age-structured multi-stock population dynamics model</a>. 

Code development is occurring, under Mac OS 10.9.5, using <a href="https://gcc.gnu.org/wiki/GFortran" target="_blank">the free, open source, GNU Fortran 95 compiler</a>. In the near future, a related open-source repository will be linked with R code to visualize and analyze the output files from the Fortran code. 

Links are provided at the top of this page to view the files in this repository, as well as to download the files therein (raw code or otherwise). Motivation and background for this project is provided below, but first a quick FAQ.

## FAQ
__Q1:__ 

What does "open source" meen?

__A1:__ 

You have the power to fix the typo in __Q1__. 

All of the files in this repository (Fortran, R, markdown, or otherwise) are available to the public, including the markdown text file underlying this webpage. Git provides us a formal process of version control for collaborating on these files. 

If you're new to Git and version control, and interested in learning how to collaborate on GitHub, <a href="https://www.youtube.com/watch?v=SCZF6I-Rc4I" target="_blank"> here's a nice video (6min 33s) with an example tututorial. </a> As a first exercise, you could access  <a href=https://github.com/John-Brandon/PBR-Tier-System/tree/master target="_blank"> the markdown text for this webpage in the README.md file </a>, following the instructions in the video. Once you've edited the typo above, submit your first pull request. There's no better way to learn than by doing. 

Regardless of whether you know how to use Git, if you see something that is not factually correct; important information is missing; find a bug in the code; have a better way of doing something -- let us know (<a href="https://guides.github.com/features/issues/" target="_blank">here's one avenue through GitHub</a>). You can also contact the maintainer of this repository, listed at the end of this page. 

Legally speaking, for this repository, open source means that the files here are licensed under the <a href="https://github.com/John-Brandon/PBR-Tier-System/blob/master/LICENSE.md" target="_blank">MIT License</a>. You can use them however you want, even for profit, but you can't hold us liable for any damages. 

__Q2__: 

Why is open source important?

__A2__: 

One of the goals of this project is to provide an opportunity to engage and grow the scientific community for simulation research on "Potential Biological Removal" (*aka* PBR -- More on PBR below). We've decided to use GitHub as an open source platform to ensure that this research is transparent and reproducible. 

We hope this project excites ideas and methodological approaches for reducing scientific uncertainty in the context of meeting the management objectives of the U.S. Marine Mammal Protection Act, or other related conservation objectives. We also offer this as a platform to test existing strategies to calculate PBR, which haven't been formally tested yet. For example, averaging abundance and mortality estimates (Wade and Demaster 1999; <a href="http://www.nmfs.noaa.gov/pr/pdfs/sars/gamms2005.pdf" target="_blank">NMFS 2005</a>) and setting recovery factors for endangered species (<a href="http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.221.9577&rep=rep1&type=pdf" target="_blank">Taylor et al. 2003</a>).


PBR is also a useful concept for scientists outside the U.S., who are looking to provide management advice on sustainable levels of human caused mortality. An open-source repository engenders international collaboration and improvements in scientific advice.  

Finally, you just never know who might be reading this and want to collaborate on a scientifically interesting and challenging problem. For all we know, the reader might be a curious computer scientist working for a tech company, with a desire to learn more about, and perhaps engage in computational applications in natural resource science. Making this project open source allows anyone with an internet connection to contribute to the code and development of algorithms.

With all of that said, let's get into the fun, starting with some background...

## The U.S. Marine Mammal Protection Act and Potential Biological Removal
The <a href="http://www.nmfs.noaa.gov/pr/laws/mmpa/text.htm" target=_blank>U.S. Marine Mammal Protection Act (MMPA)</a> contains rules to limit the number of human caused mortalities to individual stocks of marine mammals (stocks are the unit to be managed). In particular, the legislation prescribes a "Potential Biological Removal" limit to human caused mortality. Further, the MMPA mandates that this limit should be calculated as the product of three values: (1) An estimate of abundance that "provides reasonable assurance that the stock size is equal to or greater than the estimate"; (2) One half of the maximum intrinsic population growth rate; and (3) A recovery factor between 0.1 and 1.0. 

PBR can therefore be written as the equation: 
PBR = N<sub>MIN</sub> 0.50 R<sub>MAX</sub> F<sub>R</sub> 

The relevant management objective of PBR for this project, as outlined in the MMPA and quantified by the U.S. National Marine Fisheries Science Service (NMFS), is to allow stocks of marine mammals to be maintained at an abundance that is 50% or more of the carrying capacity of their environment. 

## Previous simulation research on fisheries management strategies 
The earliest examples of simulation research on fisheries management strategies were by Southward (1968) and Hilborn (1979). 

The concept later gained traction and made significant strides through the Scientific Committee of the International Whaling Commission (the SC). The SC's research and development started during the 1980s and is ongoing (<a href="http://icesjms.oxfordjournals.org/content/64/4/603.full.pdf+html" target=_blank>see Punt and Donavon, 2007</a>, for a review]). 

Generally speaking, MSE procedures that have been developed (by the SC and others since) take into account scientific uncertainties as represented by a range of 'what if' scenarios. For example, "What if abundance estimates are biased?" Would the catch control rule still meet the management objectives? Each scenario is known as a "trial". Simulation modeling is used to evaluate the performance of alternative management strategies, across a set of trials, given quantified management objectives. 

This approach for providing scientific advice to fisheries management has since become increasingly popular in recent decades, and is now also used for evaluating management systems of commercial fisheries for fin- and shell-fish. The procedures originally developed by the SC are now known in fisheries science as "Management Strategy Evaluation" (MSE). __TODO__ (REF recent review(s) e.g., Smith et al Punt et al, in fisheries)

## Management strategy evaluation 
MSE involves a few key steps, including: (i) development of an "operating model", which represents the system being managed (the 'true' state of nature) and how observable data are generated; (ii) identification of a candidate management strategy, which also includes a data collection strategy (management strategies do not exist in a vacumn re: survey effort);  (iii) developing a set of trials that span the range of key uncertainties about the true state of nature, including trials where data are biased and/or imprecise; and (iv) evaluation, using simulation, of the candidate management strategies. 

Candidate management strategies are evaluated based on their performance relative to the quantified management objectives (e.g., maintaining abundance above some level with 95% certainty). In the context of this project, the management objectives of the U.S. Marine Mammal Protection Act are the benchmark. 

Modeling is never going to be a perfect science. At best, it is a caricture of reality founded on statistics.  
<a href="http://www.fao.org/fishery/eaf-net/eaftool/eaf_tool_50" target=_blank>The Food and Agricultural Organization (FAO) of the United Nations provides a summary and critique here of the MSE simulation approach for providing scientific advice to fisheries management.</a> The FAO also notes some serious pragmactic barriers to entry re: MSE (e.g., start-up costs, access to individuals with experience in simulation modeling, etc.). 

## Previous simulation research on PBR
NMFS scientists are active members of the SC. There is cross-pollination of ideas (when things work well) between scientists from different countries working together within the SC. Initial work on simulation testing of PBR fruited from the MSE procedures being developed by the SC during the late 1980s and early 1990s. 

[Taylor (1993)](https://swfsc.noaa.gov/publications/TM/SWFSC/NOAA-TM-NMFS-SWFSC-188.PDF) used simulation to explore the impact of assigning to N<sub>MIN</sub> alternative percentiles of an abundance estimate, given survey sampling error (e.g. assigning the mean vs lower 95th confidence interval of the abundance estimate to N<sub>MIN</sub>). Taylor (1993) demonstrated that a strategy of assigning the point estimate of abundance to N<sub>MIN</sub> resulted in poor performance when biase existed in the abundance data. A strategy of assigning the lower 95th confidence interval of the abundance estimate to N<sub>MIN</sub> was relatively robust to the biases that were explored. 

Wade (1998) expanded on the simulation studies by Taylor (1993). As a first step, Wade (1998) iterated over a set of trials where estimates of abundance were unbiased, but in some cases imprecise. For these trials, Wade (1998) essentially solved for the percentile of the abundance estimate that resulted in meeting the management objective of the MMPA in 95% of the simulations across scenarios. That value was the lower 20th percentile of the abundance estimate. Wade (1998) then examined a variety of trials where different biases existed in the data. For these trials, N<sub>MIN</sub> was set to the lower 20th percentile of the abundance estimate, and the process was repeated to find the solution for the value of F<sub>R</sub> that resulted in 95% of the simulations meeting the management objectives if data are biased. The results of Wade's simulation testing form the basis for the values NMFS currently uses when calculating PBR for stocks of marine mammals. 

Below, we describe how the PBR Tier System study aims to expand on previous research. We also provide some rationale for doing so. 

### Alternative tiers of data availability <-> alternative management strategies

__TODO__: Add some general rationale for a tier system, before getting into modeling. Could link the poker example here (draft a flight simulator analogy, poker might be esoteric?): https://github.com/John-Brandon/PBR-Tier-System/blob/gh-pages/poker.md

The Tier System aspect of this research involves different levels of data availability. Previous PBR research has largely  assumed that PBR would be calculated based on a single estimate of absolute abundance. 

(although see Wade and DeMaster 1998, for an example of a PBR simulation study where multiple estimates were averaged: __TODO__ Add this reference and wrap up the remaining last bits of text)

In this case, features of the proposed tier system might include how historical abundance estimates are weighted, how trends are estimated, and whether abundance data older than 8 years are used.

__TODO__: Add image from slides showing hierarchy of tiers. 

## Operating model synopsis
The underlying population dynamics were represented using an age-, sex-, and space-structured multiple-stock model (which also allows single-stock scenarios to be evaluated). The age-structured component of the model mirrors single-stock assessment models employed by the NMFS scientists and the SC of the IWC (e.g., Punt 1999, Wade 2002). The single-stock model is expanded here to allow for two spatially overlapping stocks of the same species (Fig 1). This approach is similar to multi-stock models developed by the SC of the IWC for MSEs under more complicated scenarios, involving overlapping stocks, and uncertainties in stock identity with respect to catch limits (e.g. IWC 2012).

![Spatial Structure of the Model](https://dl.dropboxusercontent.com/u/31861309/SpatialFig.png)

*Figure 1. The ranges and overlap of the two stocks are shown over the four areas. The surveyed area comprises areas 1, 2, and 3. Estimates of abundance include all animals (both stocks) in the survey area. The ranges of the two stocks overlap in area 2. If the operating model is run as a single stock, the model collapses to the box joining areas A and B.* 

The modeling approach differs from the age-aggregated model used in previous PBR research (e.g., Taylor 1993; Wade 1998). For simple scenarios, the two approaches should result in equivalent dynamics. However, the more complex population dynamics model developed here allows a wider spectrum of trials to be evaluated (e.g., spatial overlap of stocks of the same species, and non-uniform age- and sex-specific human caused mortality).

__TODO__: Add preliminary image showing population trajectories. 

### Funding and Contributors
This project is funded by the [U.S. Western Pacific Fisheries Management Council](http://www.wpcouncil.org/about-us/). Drs. John R. Brandon, [Andr&eacute; E. Punt](http://fish.washington.edu/people/punt/index.html), Paula Moreno and Randall R. Reeves are collaborators. They are working together as [the Independent Advisory Team on marine mammal stock assessment] (http://www.usm.edu/gcrl/scemfis/marine.mammal.asssessment.team.php), with the goal of developing research projects that serve to reduce uncertainty in marine mammal stock assessment. The team is funded through [Science Center for Marine Fisheries](http://scemfis.org/aboutus.html). The center is run under the [U.S. National Science Foundation's Industry and University Cooperative Research Program.](http://www.nsf.gov/eng/iip/iucrc/program.jsp) 

### References 
Hilborn, R. 1979. Comparison of fisheries control-systems that utilize catch and effort data. *J. Fish. Res. Board. Can.* 36:1477–1489. 

IWC. 2012. International Whaling Commission. 2012. Report of the Scientific Committee. Annex E. Report of the Standing Working Group on an Aboriginal Subsistence Whaling Management Procedure. *J. Cet. Res. Manage.* (Suppl.) 13:130-53. 

[NMFS. 2005. Revisions to Guidelines for Assessing Marine Mammal Stocks. 24 pp. Available at: http://www.nmfs.noaa.gov/pr/pdfs/sars/gamms2005.pdf](http://www.nmfs.noaa.gov/pr/pdfs/sars/gamms2005.pdf)

Punt, A. E. 1999. A full description of the standard Baleen II model and some variants thereof. *J. Cet. Res. Manage.* 1:267-276

Southward, G. M. 1968. A simulation of management strategies in the Pacific halibut fishery. Technical Report 47, International Pacific Halibut Committee, Seattle, WA.

[Taylor, B.L. 1993. “Best” abundance estimates and best management: Why they are not the same. U.S. Department of Commerce, NOAA Technical Memorandum NMFS-SWFSC-188. 20 pp.](https://swfsc.noaa.gov/publications/TM/SWFSC/NOAA-TM-NMFS-SWFSC-188.PDF)

[Taylor, B.L., Scott, M., Heyning, J. and Barlow, J. 2003. Suggested guidelines for recovery factors for endangered
marine mammals. U.S. Dep. Commer., NOAA Tech. Memo. NOAA-TM-NMFS-SWFSC-354, 6 pp.](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.221.9577&rep=rep1&type=pdf)

Wade, P.R. 1998. Calculating limits to the allowable human-caused mortality of cetaceans and pinnipeds. *Mar. Mamm. Sci.* 14: 1-37.

Wade, P.R. 2002. A Bayesian stock assessment of the eastern Pacific gray whale using abundance and harvest data from 1967-1996. *J. Cet. Res. Manage.* 4: 85-98. 

Wade, P.R. and DeMaster, D. P. 1999. Determining the optimimum interval for abundance surveys. Marine Mammal Survey and Assessment Methods, Garner *et al. (eds)*. Balkema, Rotterdam

---------------------------

*last updated: 08 Aug 2015*

*contact: jbrandon@gmail.com*
