# Data files and code accompanying Shilts _et al._, 2022

This repository contains all of the data files and code needed to reproduce the analyses published in the paper <br>
__A physical wiring diagram of the human immune system__ <br>
Jarrod Shilts _et al._, 2022
<br>
<br>
When reusing any of these datasets or scripts, please cite our paper : [URL](https://www.nature.com/articles/s41586)
<br>

## What is included here

Each folder contains both the raw and processed data files that were used in this study, along with the custom scripts that run the paper's calculations and plot figures.
For most of the data sets that are most likely to be useful to re-use by others, we have already provided these as Supplementary Tables, which can be accessed at the [journal's website](https://www.nature.com/articles/s41586#Extended). Therefore these files are provided predominantly for more advanced projects that aim to build off our study, as well as to ensure reproduciblity and transparency. 
<br><br><br>
Given the diverse types of analyses included in this larger study, the repository is divided into folders that each cover a different section of the paper: <br><br>
__organized_code_screen_processing__ : raw and processed data from the comprehensive protien-protein interaction screening done using the SAVEXIS technique to discover receptor binding partners. The code covers the processing of screen data and benchmarking. This mostly covers Figure 1 of the paper. <br><br>
__organized_code_affinity__ : datasets of surface interactions integrated with proteomics expression. The code analyzes these datasets to test hypotheses about interaction abundance and binding kinetics. This mostly covers Figure 2 of the paper. <br><br>
__organized_code_modelling__ : parameter values for the mathematical model of cell-cell association, along with the associated code for evaluating the model. This covers parts of Figures 2 and Figure 4 of the paper. <br><br>
__organized_code_scRNA_website__ : integrated datasets of single-cell RNAseq and surface protein interaction matrices. The code provides a range of interactive plotting funcitons.  This covers the first part of Figure 3 of the paper.  <br>
This particular code base can be explored more conveninently for day-to-day use by accessing our web tool :https://www.sanger.ac.uk/tool/immune-interaction/immune-interaction/  <br><br>
__organized_code_scRNA_analysis__ : an extension of the folder above that contains additional code covering the middle sections of Figure 2. <br><br>
__organized_code_spatialRNA__ : spatial transcriptomics data from human lymph node provided by 10x Genomics. The code calculates the spatial colocalization of gene pairs. This covers the last portion of Figure 3. <br><br>
__organized_pharm_analysis__ : high-content microscopy data of isolated human immune cells stimulated with recombinant surface proteins. The code processes the data and plots the phenotypes found from the experiment. This covers Figure 4 of the paper. <br><br>


### How to run the code

Code is written in either the Python or R programming languages. A summary of the versions and package depencies  needed to run the code are summarized in the __versions.txt__ file. For conveninence, the R packages are also provided as a __renv.lock__ file.


