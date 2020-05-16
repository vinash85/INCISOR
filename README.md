###
The datasets to run code is available at  http://www.umiacs.umd.edu/~vinash85/public/incisor.tar.gz 
The working directory (working dir) contains two subdirectories, src and data, which respectively contain the code and the input data set.
The code needs to be run from working directory. 
###  
The  src subdirectory contains two versions of the INCISOR code:
1. INCISOR.DU.R: The Full code INCISOR mining all 500 million gene-gene pairs to identify cancer SR interactions. The code can be run by simply sourcing the INCISOR.DU.R.  Once the code run successfully the final output, i.e the list of SR interaction identified by INCISOR, will be stored in the variable "positive.sr.du.interactions". This variable also stores the aggregate score/significance obtained at each of the four screening steps of INCISOR.  All the input data set required are automatically uploaded by the R code from “data” subdirectory . 

2. INCISOR.DU.top500.R: This is a version INCISOR that runs on top 500 SR interaction gene pairs that were identified via applying the full INCISOR code. In particular, it runs all four screening steps (in-vitro, sof, clinical and phylogenetic screens). Please note that FDR level estimated on this subset will be obviously  different from the full INCISOR that processes all 500 million gene pairs given the much smaller search space. The final output, positive SR interactions and corresponding score/signifcance will be stored in a variable named "positive.sr.du.interactions".  All the input data set required are automatically uploaded by the R code from “data” subdirectory. 
As a negative example of SR interactions, we chose 500 gene pairs that show shRNA evidence of interactions but were not identified by INCISOR having SR interactions. Set negative control flag to TRUE (use.negative.sr.interaction=TRUE) in INCISOR.DU.top500.R to activate this mode.  The result of each screens of INCISOR for negative examples are given in variable "sr.gene.all.screen.stats"   

The src subdirectory also contains:
source.incisor.R: This file contain main helper functions of INCISOR, that are needed to run it.
HyperGeometricTest.pair.cpp: C-code that is sourced by (using Rcpp) both the above versions of INCISOR. 


### README ###

INCISOR need following R library installed:
data.table
Rcpp
RcppArmadillo
Parallel
survival
lsr
doMC
foreach

#### In some environment to run the c-code following commands are required to set the variables 
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
library(RcppArmadillo)




