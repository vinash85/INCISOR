##### INCISOR source code #####

############# load libraries
library(data.table)
library(Rcpp)
library(parallel)
library(survival)

require(doMC)
require(foreach)
registerDoMC(cores = 64)

####### load TCGA data ##########
load("TCGA.RData")

####### The mRNA expression (RNAseqV2) and patients' clinical characteristics were downloaded from TCGA Data Portal,
### and copy number data (Gistic values) was downloaded from Broad Institute's Firehose (https://gdac.broadinstitute.org/)
### on June 27, 2015. The TCGA data includes the following data fields:
### genes         19001 -none-     character: protein coding genes' symbols
### entrez        19001 -none-     numeric  : protein coding genes' entrez ID
### samples        8749 -none-     character: TCGA sample barcodes
### types          8749 -none-     character: cancer types 
### mRNA      166239749 -none-     numeric  : a matrix of gene expression (RNAseq) with genes in the row and samples in the column
### scna      166239749 -none-     numeric  : a matrix of SCNA (Gistic values) with genes in the row and samples in the column
### mRNAq2    166239749 -none-     numeric  : a matrix of mRNA expression equally divided into 3 groups (low:0, middle:1, high:2) using 1/3, 2/3-quantile across the samples in each cancer type
### scnaq2    166239749 -none-     numeric  : a matrix of SCNA equally divided into 3 groups (low:0, middle:1, high:2) using 1/3, 2/3-quantile across the samples in each cancer type
### sex            8749 -none-     character: sex of patients
### age            8749 -none-     numeric  : age of patients 
### race           8749 -none-     character: ethnicity of patients
### survival      17498 -none-     numeric  : two column matrix with first column as survival time (in days) and second column as cencering (death=0, alive=1)
### stage          8749 -none-     character: stage of the tumor 
### mRNA.norm 166239749 -none-     numeric  : quantile-normalized of mRNA expression values
### scna.norm 166239749 -none-     numeric  : quantile-normalized of SCNA values
### surv.dt           2 data.table list		: survival matrix for Cox regression (first column: survival time (in days), second column: censoring (death=1, alive=0))

### all the data inside the pancancer data.structure. 
pancancer = prob 

source("source.incisor.R")

############# store variables for following analysis
genes=pancancer$genes
numGenes=length(genes)
numSamples=length(pancancer$samples)
cancerType = pancancer$types;typeInx = as.numeric(as.factor(cancerType))
typeNum= length(unique(typeInx ));levels(typeInx) = seq(typeNum)-1
mRNA.norm = pancancer$mRNA.norm
scna.norm = pancancer$scna.norm
surv.dt = pancancer$surv.dt
surv.type=as.numeric(as.factor(pancancer$types))
age=qnorm.array(pancancer$age)
race=pancancer$race
sex=pancancer$sex

##########################################################################################
################### Collecting P-values from the 4 INCISOR screenings #######################
##########################################################################################

##########################################################################################
###### Step 1. SoF screening ############################################################
##########################################################################################
sourceCpp("HypergeometricTest.cpp")
numGenes = nrow(scna)
numSamples = ncol(scna)

xs=seq(numGenes)

sourceCpp("HyperGeometricTest.pair.cpp")

########## molecular scrneeing (hyperGeometric  test) for SCNA:
########## outbed[,1:2]: gene indices, outbed[,3:4]: the depletion of non.rescued state. oubed significance of gene A-gene B low depletion
########## negative p-values -> depletion, positive p-values -> enrichment
#### given a gene activity in (0: inactive, 1 : normal and 2: high), there are 9 possible functional activity state of gene pair A-B. 

outbed = foreach(x = xs, .inorder=T, .combine=rbind) %dopar%{
	scnacurr = scnaq[x,] 
	scnacurr[is.na(scnacurr)] =1
	out = hypergeometricTest(scnacurr,scnaq)## conduct hypergeometric test for all possible 9 functional states of a gene pair. 
	out[,1] = ifelse(out[,1] < 0, out[,1], 1-out[,1]) ## rescued state enrichment p.value converted to depletion.  
	out[,2] = ifelse(out[,2] < 0, out[,2], 1-out[,2]) ## rescued state enrichment p.value converted to depletion.  
	out2=cbind(x,seq(numGenes),out)
	print(x)
	return(out2)
}
all.pairs.tested = outbed[,1:2]
#### note in outbed NA is represented as p.value < 0 
molecular.scna.fdr = sapply(1:3, function(tt) {
	p = ifelse(outbed[,tt+2] < 0 , 1, outbed[,tt+2])
	p.adjust(p, method="fdr")
})
##for mRNA
outbed = foreach(x = xs, .inorder=T, .combine=rbind) %dopar%{
	mRNAcurr = mRNAq[x,] 
	mRNAcurr[is.na(mRNAcurr)] =1
	out = hypergeometricTest(mRNAcurr,mRNAq)
	out[,1] = ifelse(out[,1] < 0, out[,1], 1-out[,1]) ## rescued state enrichment p.value converted to depletion.  
	out[,2] = ifelse(out[,2] < 0, out[,2], 1-out[,2]) ## rescued state enrichment p.value converted to depletion.  
	out2=cbind(x,seq(numGenes),out)
	print(x)
	return(out2)
}
#### note in outbed NA is represented as p.value < 0 

molecular.mRNA.fdr = sapply(1:3, function(tt) {
	p = ifelse(outbed[,tt+2] < 0 , 1, outbed[,tt+2])
	p.adjust(p, method="fdr")
})

molecular.fdr.thr = 0.05  # 10% 
molecular.rescued.screen = which( 
	molecular.scna.fdr[,3] < molecular.fdr.thr &
	molecular.mRNA.fdr[,3] < molecular.fdr.thr 
	)
molecular.non.rescued.screen = which( 
	(molecular.scna.fdr[,1] < molecular.fdr.thr |
	molecular.scna.fdr[,2] < molecular.fdr.thr )&
	(molecular.mRNA.fdr[,1] < molecular.fdr.thr |
	molecular.mRNA.fdr[,2] < molecular.fdr.thr )
	)

sof.screen =  intersect(molecular.rescued.screen,molecular.non.rescued.screen )


sr.sof = all.pairs.tested[sof.screen,] ### putative SR gene pairs passing SOF screen 


##########################################################################################
###### Step 2. Clinical screening ########################################################
##########################################################################################
### pre-processing the molecular and survival data for cox regression 

sr.curr = sr.sof
library(parallel)
mRNA.norm = mclapply(1:numSamples, function(tt) qnorm.array(pancancer$mRNA[,tt]),mc.cores=64)
scna.norm = mclapply(1:numSamples, function(tt) qnorm.array(pancancer$scna[,tt]),mc.cores=64)
mRNA.norm = t(do.call(rbind, mRNA.norm))
scna.norm = t(do.call(rbind, scna.norm))

mRNA.norm = mclapply(1:numGenes, function(tt) qnorm.array(mRNA.norm[tt,]),mc.cores=64)
scna.norm = mclapply(1:numGenes, function(tt) qnorm.array(scna.norm[tt,]),mc.cores=64)
mRNA.norm = do.call(rbind, mRNA.norm)
scna.norm = do.call(rbind, scna.norm)
scnaq2 = pancancer$scnaq2
mRNAq2 = pancancer$mRNAq2


pancancer$age = ifelse(is.na(pancancer$age), mean(pancancer$age,na.rm=T), pancancer$age)
pancancer$race = ifelse(is.na(pancancer$race),"unknown", pancancer$race)
pancancer$sex = ifelse(is.na(pancancer$sex),"FEMALE", pancancer$sex)
surv.all = data.table(pancancer$survival, type= pancancer$types, age = qnorm.array(pancancer$age), sex = pancancer$sex, race = pancancer$race)
setnames(surv.all, 1:2, c("time","status"))
surv.all$status = ifelse(surv.all$status==1,0,1)
pancancer$surv.strata = surv.all

##########################################################################################
### cox regression using scna 
##########################################################################################
cox.scna = mclapply(1:nrow(sr.curr), function(tt) cox.pair.du.strata.controlled(sr.curr[tt,], prob =pancancer, use.mRNA=F),mc.cores=64))
cox.du.scna.curr = do.call(rbind, cox.scna)
cox.du.scna.curr=cbind(sr.curr,cox.du.scna.curr)

##########################################################################################
### cox regression using mRNA 
##########################################################################################
cox.mRNA = mclapply(1:nrow(sr.curr), function(tt) cox.pair.du.strata.controlled(sr.curr[tt,], prob =pancancer, use.mRNA=T),mc.cores=64))
cox.du.mRNA.curr = do.call(rbind, cox.mRNA)
cox.du.mRNA.curr=cbind(sr.curr,cox.du.mRNA.curr)

clinical.beta =(cbind( cox.du.mRNA[,c(3,9)], cox.du.scna[,c(3,9)]))
class(clinical.beta) = "numeric"
clinical.beta[,c(2,4)] = - clinical.beta[,c(2,4)] ### beta negative are desirable from here on 
clinical.p = cbind( cox.du.mRNA[,c(7,13)],cox.du.scna[,c(7,13)] )
class(clinical.p) = "numeric"

clinical.fdr = sapply(1:ncol(clinical.beta), function(tt) {
	p = ifelse(clinical.beta[,tt] < 0 , clinical.p[,tt], 1 )
	p.adjust(p, method="fdr")
})


clincial.fdr.thr = .05
clinical.screen = which(rowSums(clinical.empirical.fdr < clinicalfdr.thr ) >= 4 ) 

sr.clincial = sr.curr[clinical.screen,]


##########################################################################################
###### Step 3. shRNA screening ######################################################
##########################################################################################
### We collected the following 4 large-scale shRNA essentiality screening datasets for step 3, 
### and run the screening for all possible gene pairs
load("achilles.old.RData") # Cheung et al. PNAS (2011).
load("achilles.new.RData") # Cowley et al. Sci. Data. (2014).
load("marcotte.old.RData") # Marcotte et al. Cancer Discov (2012).
load("marcotte.new.RData") # Marcotte et al. Cell (2016).

### this table shows the data fields included in the first dataset (Cheung et al.) as an example,
### the remaining 3 datasets have a similar data structure.
###          Length  Class      Mode     
###genes       19001 -none-     character : gene symbols of protein coding genes
###entrez      19001 -none-     numeric   : entrez IDs of protein coding genes
###samples       102 -none-     character : sample IDs (as coded in the dataset)
###celllines     102 -none-     character : names of cell lines
###types         102 -none-     character : cancer types of cell lines
###mRNA      1938102 -none-     numeric   : matrix of mRNA expression with genes in rows, cell lines in columns
###mRNAq     1938102 -none-     numeric   : matrix of mRNA expression equally divided into 3 groups (low:0, middle:1, high:2) using 1/3, 2/3-quantile across cell lines
###scna      1938102 -none-     numeric   : matrix of SCNA expression with genes in rows, cell lines in columns
###scnaq     1938102 -none-     numeric   : matrix of SCNA expression equally divided into 3 groups (low:0, middle:1, high:2) using 1/3, 2/3-quantile across cell lines
###mat           102 data.frame list      : matrix of shRNA essentiality scores with genes in rows, cell lines in columns (the lower the more essential)



### the following code shows how to obtain the shRNA screeing significance of all gene pairs
### for the first dataset (Cheung et al.) as an example
#########################################################################################################
### one can repeat the process for the remaining 3 datasets to obtain p.ach.new, p.mar.old, p.mar.new ###
#########################################################################################################
sr.curr = sr.clincial 
ess = prob 
sr.ess = cbind( match(pancancer$genes[sr.curr[,1]],ess$genes), match(pancancer$genes[sr.curr[,2]],ess$genes))
require(doMC)
require(foreach)
registerDoMC(cores = 64)
outbed = foreach(x = 1:nrow(sr.ess), .inorder=T) %dopar% pair.ess.wil(sr.ess[x,], ess = ess)
shrna.wil.mat = do.call(rbind, outbed)
shrna.wil.mat.fdr = apply(shrna.wil.mat, 2, p.adjust, method="fdr")

shrna.fdr.thr = 0.1 

rescuer.screen = rowSums(cbind(shrna.wil.mat.fdr[,c(1,3)], shrna.wil.bc.mat[,c(1,3)]) < shrna.fdr.thr , na.rm=T) >=1
vulnerable.screen = rowSums(cbind(shrna.wil.mat.fdr[,c(2,4)], shrna.wil.bc.mat[,c(2,4)]) < shrna.fdr.thr , na.rm=T) >=1
shrna1.screen = (rescuer.screen & vulnerable.screen) 

### similar screen is obtained for all the other 3 shRNA dataset 


sr.shrna = sr.curr[shrna1.screen|shrna2.screen|shrna3.screen|shrna4.screen, ] 


##########################################################################################
###### Step 4. Phylogenetic screening ####################################################
##########################################################################################
load("yuval.phylogenetic.profile.RData")
### The phylogenetic profile is downloaded from Yuval Tabach et al. Mol Syst Biol. (2013), Supplementary Table 1
load("cluster.weight.RData")
### the feature weights are determined based on the phylogenetic tree (Ensembl database: http://useast.ensembl.org/index.html)

### function to identify phylogenetic distance of a pair of genes
sr.curr  =sr.shrna

sr.gene.all = data.table(rescuer = pancancer$genes[sr.curr[,1]], vulnerable = pancancer$genes[sr.curr[,2]])
phylo.score = phylo.profile(sr.gene.all )


### calculate the phylogenetic distance of random gene pairs for false discovery correction
sr.rand=cbind(prob$genes[sample(numGenes,100000,replace=T)],prob$genes[sample(numGenes,100000,replace=T)])
na.inx=which(sr.rand[,1]!=sr.rand[,2])
sr.rand=sr.rand[na.inx,]
pp.rand=phylo.profile(sr.rand)

########## Select pairs that pass Step 4: Phylogenetic screening
phylogenetic.screen =which( phylo.score <quantile(as.numeric(pp.rand),0.1,na.rm=T)))
sr.phylogenetic = sr.curr[phylogenetic.screen, ]

##########################################################################################
########## Final SR pairs that pass all 4 screenings #####################################
##########################################################################################

sr.du.final= sr.phylogenetic











