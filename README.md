# INCISOR

Identifying synethetic rescue interactions in humans. 

Check out the manuscript in Molecular System Biology: 
- [msb website](https://www.embopress.org/doi/full/10.15252/msb.20188323)

# System requirements 

INCISOR has been tested on R versions >= 3.2. INCISOR has been tested on Linux and OS X.

# Installation  

First download the dataset from https://www.dropbox.com/s/uw22lmqvmkkwezn/incisor.tar.gz?dl=0
The working directory (working dir) contains data directory containing the input data set.

Open R from INCISOR working directory 

```
gunzip -c incisor.tar.gz|tar -xvf
mv incisor/data INCISOR/
cd INCISOR 
R
```
To run INCISOR, install dependencies
```r
install.packages(c("Rcpp", "RcppArmadillo", "survival", "lsr", "doMC", "parallel", "data.table", "foreach"))
```
# Extra settings
In some environment to run the c-code following commands are required to set the variables : 

```r
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
library(RcppArmadillo)
```

## Quick start 

```r
source("src/INCISOR.DU.top500.R")
```
 "INCISOR.DU.top500.R" identifies top 500 SR interaction gene pairs. It runs all four screening steps (in-vitro, sof, clinical and phylogenetic screens). FDR level estimated on this subset will be obviously  different from the full INCISOR. The final output is "positive.sr.du.interactions".  

## Screen all 500 million gene pair 
Note the code takes 3 days in a 32 core machines
```r
source("src/INCISOR.DU.R")
```
Results are available in "positive.sr.du.interactions"






