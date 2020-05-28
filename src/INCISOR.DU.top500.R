##### INCISOR source code #####

use.negative.sr.interaction=FALSE 
### set this flag to TRUE to include 500 negative example of SR interaction. 
# use.negative.sr.interaction = TRUE

############# set ENV and load librariesa
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(survival)
library(lsr)
require(doMC)
require(foreach)
registerDoMC(cores = 64)

####### load TCGA data ##########
load("data/TCGA.RData")

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

### all data inside the pancancer data.structure. 
pancancer = prob 

source("src/source.incisor.R")

############# store variables for following analysis
genes=pancancer$genes
numGenes=length(genes)
numSamples=length(pancancer$samples)
cancerType = pancancer$types;typeInx = as.numeric(as.factor(cancerType))
typeNum= length(unique(typeInx ));levels(typeInx) = seq(typeNum)-1

##########################################################################################
################### Collecting P-values from the 4 INCISOR screenings #######################
##########################################################################################
### top 500 SR interactions 
 
if(use.negative.sr.interaction){
	load("data/sr.neg.500.RData")
	sr.all = sr.neg.500
}else{
	load("data/sr.500.RData")
	sr.all = sr.500
}
# sr.all = rbind(rep(seq(numGenes), each=numGenes), rep(seq(numGenes), numGenes))
# sr.all = sr.all[sr.all[,1]!=sr.all[,2],]

##########################################################################################
###### Step 2. invitro  screening ######################################################
##########################################################################################
### We collected the following 2 large-scale shRNA (achilles and marcotte) and a drug screen essentiality screening datasets for step 1, 
### and run the screening for all possible gene pairs






##### drug response screen (Iorio et. al 2016, Cell)

iorio = local({load("data/iorio.RData");environment()}); iorio =iorio$prob;
############################################
### dataset structure is same as of essentiality screen dataset  described above with exceptions of matrix named "mat" now contain IC50 values 
###mat           102 data.frame list      : matrix of ic50  of drugs in rows, cell lines in columns (the lower the more essential)
#### it also contains additional 
#               Length   Class  Mode     
# celllines          628 -none- character
# drugs              221 -none- character: list of drugs screened 
# drugTargetMap      862 -none- character: a data.frame with mapping between drug and gene targets, with first column as drug and second column as it target. 
############################################

###identify best reversal drug - gene pairs. 
iorio$mRNAq = t(apply(iorio$mRNA, 1,  function(tt) {
	probs  = quantile(tt, c(.1),na.rm=T)
	(tt > probs) + 0
}))
iorio$scnaq = t(apply(iorio$scna, 1,  function(tt) {
	probs  = quantile(tt, c(.1),na.rm=T)
	(tt > probs) + 0
}))


iorio$mat = apply(iorio$mat, 2, as.numeric)

library(effsize)
require(doMC)
require(foreach)
registerDoMC(cores = 64)
outbed = foreach(x = seq(length(iorio$genes)), .inorder=T) %dopar% find.interacting.sig(x,iorio=iorio, probs=0.1, alternative="less")
iorio.sig = do.call(cbind, outbed) ### s 
outbed.eff = foreach(x = seq(length(iorio$genes)), .inorder=T) %dopar% find.interacting.eff(x,iorio=iorio, probs=0.1, alternative = "less")
iorio.eff = do.call(cbind, outbed.eff) ###  

fdr.thr.curr = 0.05 
iorio.sig.fdr = t(apply(iorio.sig, 1, function(tt) p.adjust(tt, method="fdr")))

iorio.sig.pairs.large = which(iorio.sig.fdr < fdr.thr.curr, arr.ind = T)
iorio.temp = mclapply(1:nrow(iorio.sig.pairs.large), function(tt){
	drug.curr = iorio$drugs[iorio.sig.pairs.large[tt,1]]
	gene.curr = iorio$drugTargetMap[iorio$drugTargetMap[,1]==drug.curr,2]
	target.curr = which(pancancer$genes %in% gene.curr)
	out =NULL
	if(length(gene.curr) > 0) out = cbind(target.curr,iorio.sig.pairs.large[tt,2] )
		out 
}, mc.cores =64)
iorio.sr.pos = do.call(rbind, iorio.temp)

iorio.sig.pairs.gene = data.table(
	drug = iorio$drugs[iorio.sig.pairs.large[,1]],
	gene = pancancer$genes[iorio.sig.pairs.large[,2]],
	score = iorio.eff[which(iorio.sig.fdr < fdr.thr.curr )]
	)


  temp1 = mclapply(seq(nrow(iorio.sig.pairs.gene)), function(tt) find.gene.pair.iorio(iorio.sig.pairs.gene[tt]), mc.cores=64)

  	iorio.synergistic = do.call(rbind, temp1)
  	iorio.synergistic[,label:= paste(gene, target, sep=":")]
  	iorio.syn.score =iorio.synergistic[, list(score=max(score)), by=label]




#  iorio scna 

outbed = foreach(x = seq(length(iorio$genes)), .inorder=T) %dopar% find.interacting.sig.scna(x,iorio=iorio, alternative="less", low.thr=-1, high.thr =0)
iorio.sig = do.call(cbind, outbed) ### s 
outbed.eff = foreach(x = seq(length(iorio$genes)), .inorder=T) %dopar% find.interacting.eff.scna(x,iorio=iorio, alternative = "less", low.thr = -1, high.thr =  0)
iorio.eff = do.call(cbind, outbed.eff) ###  
iorio.sig.fdr = p.adjust(iorio.sig, method="fdr")
# fdr.thr.curr = 0.05

iorio.sig.pairs.large.scna = which(iorio.sig < fdr.thr.curr & iorio.eff >= .8,  arr.ind = T)

iorio.temp = mclapply(1:nrow(iorio.sig.pairs.large.scna), function(tt){
  drug.curr = iorio$drugs[iorio.sig.pairs.large.scna[tt,1]]
  gene.curr = iorio$drugTargetMap[iorio$drugTargetMap[,1]==drug.curr,2]
  target.curr = which(pancancer$genes %in% gene.curr)
  out =NULL
  if(length(gene.curr) > 0) out = cbind(target.curr,iorio.sig.pairs.large.scna[tt,2] )
    out 
}, mc.cores =64)
iorio.sr.pos = do.call(rbind, iorio.temp)

iorio.sig.pairs.gene.scna = data.table(
  drug = iorio$drugs[iorio.sig.pairs.large.scna[,1]],
  gene = pancancer$genes[iorio.sig.pairs.large.scna[,2]],
  score = iorio.eff[which(iorio.sig< fdr.thr.curr & iorio.eff >= .8)]
  )

temp1 = mclapply(seq(nrow(iorio.sig.pairs.gene.scna)), function(tt) find.gene.pair.iorio(iorio.sig.pairs.gene.scna[tt]), mc.cores=64)
iorio.synergistic.scna = do.call(rbind, temp1)
iorio.synergistic.scna[,label:= paste(gene, target, sep=":")]
iorio.syn.score.scna =iorio.synergistic.scna[, list(score=max(score)), by=label]

# drug.label.selected = intersect(iorio.syn.score$label,iorio.syn.score.scna$label)
drug.label.selected1 = unique(c(iorio.syn.score$label,iorio.syn.score.scna$label))

 drug.label.selected2  = sapply(drug.label.selected1, function(tt) {
	aa = unlist(strsplit(tt, split=":"))
	paste(aa[2], aa[1], sep=":")
})
  drug.label.selected = c( drug.label.selected1,  drug.label.selected2)
sr.gene.all = data.table(rescuer = pancancer$genes[sr.all[,1]], vulnerable = pancancer$genes[sr.all[,2]])

iorio.screen = paste(sr.gene.all$rescuer, sr.gene.all$vulnerable,sep=":") %in%  drug.label.selected 

temp1 = cbind( iorio.syn.score.scna$score[match(paste(sr.gene.all$rescuer, sr.gene.all$vulnerable,sep=":"), iorio.syn.score.scna$label)], 
iorio.syn.score.scna$score[match(paste(sr.gene.all$vulnerable, sr.gene.all$rescuer,sep=":"), iorio.syn.score.scna$label)], 
iorio.syn.score$score[match(paste(sr.gene.all$rescuer, sr.gene.all$vulnerable,sep=":"), iorio.syn.score$label)],
iorio.syn.score$score[match(paste(sr.gene.all$vulnerable, sr.gene.all$rescuer,sep=":"), iorio.syn.score$label)])
temp2= apply(temp1, 1, max,na.rm=T)
iorio.cohen.score = ifelse(  temp2==-Inf, NA, temp2)  

#shRNA screen 

# load("achilles.old.RData") # Cheung et al. PNAS (2011).
# load("achilles.new.RData") # Cowley et al. Sci. Data. (2014).
# load("marcotte.old.RData") # Marcotte et al. Cancer Discov (2012).
# load("marcotte.new.RData") # Marcotte et al. Cell (2016).

shRNA.screen.marcotte.new =  shRNA.screen.dataset(sr.all, shRNA.file="data/marcotte.new.RData") # Marcotte et al. Cell (2016).

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



### similar screen is obtained for all the other 3 shRNA dataset 
### the following code shows how to obtain the shRNA screeing significance of all gene pairs
### for the first dataset (Cheung et al.) as an example

shRNA.screen.marcotte.old =  shRNA.screen.dataset(sr.all, shRNA.file="data/marcotte.old.RData") # Marcotte et al. Cancer Discov (2012).
shRNA.screen.achilles.old =  shRNA.screen.dataset(sr.all, shRNA.file="data/achilles.old.RData") # Cheung et al. PNAS (2011).
shRNA.screen.achilles.new =  shRNA.screen.dataset(sr.all, shRNA.file="data/achilles.new.RData") # Cowley et al. Sci. Data. (2014).

# sum(shRNA.screen.marcotte.new < 0.05 |  shRNA.screen.marcotte.old< 0.05 | shRNA.screen.achilles.old <0.05 |  shRNA.screen.achilles.new < 0.05)


invitro.screen = iorio.screen|shRNA.screen.marcotte.new < 0.1 |  shRNA.screen.marcotte.old< 0.1 | shRNA.screen.achilles.old <0.1 |  shRNA.screen.achilles.new < 0.1 ### FDR set to 10% 

temp1 = cbind(shRNA.screen.marcotte.new,shRNA.screen.marcotte.old, shRNA.screen.achilles.old, shRNA.screen.achilles.new )
temp2= apply(temp1, 1, min,na.rm=T)
shRNA.significance.fdr = ifelse(  temp2==-Inf, NA, temp2) 

sr.curr=sr.all
sr.invitro = sr.curr[invitro.screen, ] 


##########################################################################################
###### Step 2. SoF screening ############################################################
##########################################################################################
sr.curr = sr.all
library(Rcpp)
sourceCpp("src/HyperGeometricTest.pair.cpp",rebuild=T)

########## molecular scrneeing (hyperGeometric  test) for SCNA:
########## outbed[,1:2]: gene indices, outbed[,3:4]: the depletion of non.rescued state. oubed significance of gene A-gene B low depletion
########## negative p-values -> depletion, positive p-values -> enrichment
#### given a gene activity in (0: inactive, 1 : normal and 2: high), there are 9 possible functional activity state of gene pair A-B. 

## molecular 
pval.scna1 = hypergeometricTestPair(scnaq= pancancer$scnaq2, pairs=sr.curr)
pval.scna.up = hypergeometricTestPair(scnaq= pancancer$scnaq2, pairs=sr.curr,lowerTail=0)
pval.scna = cbind(sr.curr[,1:2], pval.scna1, pval.scna.up)

pval.mRNA1 = hypergeometricTestPair(scnaq= pancancer$mRNAq2, pairs=sr.curr)
pval.mRNA.up = hypergeometricTestPair(scnaq= pancancer$mRNAq2, pairs=sr.curr,lowerTail=0)
pval.mRNA = cbind(sr.curr[,1:2], pval.mRNA1, pval.mRNA.up)


sof.fdr.scna = sof.pattern(pval.scna)
sof.fdr.mRNA = sof.pattern(pval.mRNA)
sof.aggreage.fdr = ifelse(sof.fdr.scna < sof.fdr.mRNA, sof.fdr.scna, sof.fdr.mRNA)
sof.screen = sof.fdr.scna < 0.05 & sof.fdr.mRNA < 0.05


##########################################################################################
###### Step 3. Clinical screening ########################################################
##########################################################################################
### pre-processing the molecular and survival data for cox regression 
sr.curr = sr.all

library(parallel)
numGenes = length(pancancer$genes)
numSamples = length(pancancer$sample)
mRNA.norm = mclapply(1:numSamples, function(tt) qnorm.array(pancancer$mRNA[,tt]),mc.cores=64)
scna.norm = mclapply(1:numSamples, function(tt) qnorm.array(pancancer$scna[,tt]),mc.cores=64)
mRNA.norm = t(do.call(rbind, mRNA.norm))
scna.norm = t(do.call(rbind, scna.norm))

mRNA.norm = mclapply(1:numGenes, function(tt) qnorm.array(mRNA.norm[tt,]),mc.cores=64)
scna.norm = mclapply(1:numGenes, function(tt) qnorm.array(scna.norm[tt,]),mc.cores=64)
mRNA.norm = do.call(rbind, mRNA.norm)
scna.norm = do.call(rbind, scna.norm)

pancancer$mRNA.norm = mRNA.norm 
pancancer$scna.norm = scna.norm
##imputation of missing confounders
pancancer$age = ifelse(is.na(pancancer$age), mean(pancancer$age,na.rm=T), pancancer$age)
pancancer$race = ifelse(is.na(pancancer$race),"unknown", pancancer$race)
pancancer$sex = ifelse(is.na(pancancer$sex),"FEMALE", pancancer$sex)
gii =  calculate.genomic.instability(pancancer$scna)
surv.all = data.table(pancancer$survival, types= pancancer$types, age = qnorm.array(pancancer$age), sex = pancancer$sex, race = pancancer$race, gii = gii)
setnames(surv.all, 1:2, c("time","status"))
surv.all$status = ifelse(surv.all$status==1,0,1)
pancancer$surv.strata = surv.all
pancancer$types = pancancer$type 



##########################################################################################
### cox regression using scna 
##########################################################################################
cox.scna = mclapply(1:nrow(sr.curr), function(tt) cox.pair.du.strata.controlled(sr.curr[tt,], prob =pancancer, use.mRNA=F),mc.cores=32)
cox.du.scna.curr = do.call(rbind, cox.scna)
cox.du.scna=cbind(sr.curr,cox.du.scna.curr)

##########################################################################################
### cox regression using mRNA 
##########################################################################################
cox.mRNA = mclapply(1:nrow(sr.curr), function(tt) cox.pair.du.strata.controlled(sr.curr[tt,], prob =pancancer, use.mRNA=T),mc.cores=32)
cox.du.mRNA.curr = do.call(rbind, cox.mRNA)
cox.du.mRNA=cbind(sr.curr,cox.du.mRNA.curr)


clinical.beta =(cbind( cox.du.mRNA[,c(3,9)], cox.du.scna[,c(3,9)]))
class(clinical.beta) = "numeric"
clinical.beta[,c(2,4)] = - clinical.beta[,c(2,4)] ### beta negative are desirable from here on 
clinical.p = cbind( cox.du.mRNA[,c(7,13)],cox.du.scna[,c(7,13)] )
class(clinical.p) = "numeric"

clinical.fdr = sapply(1:ncol(clinical.beta), function(tt) {
	p = ifelse(clinical.beta[,tt] < 0 , clinical.p[,tt], 1 )
	p.adjust(p, method="fdr")
})

clinical.fdr.thr = .05
clinical.screen = rowSums(clinical.fdr < clinical.fdr.thr ) >= 1 
clinical.signficance.fdr = apply(clinical.fdr, 1, min)
sr.clinical = sr.curr[which(clinical.screen),]



##########################################################################################
###### Step 4. Phylogenetic screening ####################################################
##########################################################################################
load("data/yuval.phylogenetic.profile.RData")
### The phylogenetic profile is downloaded from Yuval Tabach et al. Mol Syst Biol. (2013), Supplementary Table 1
load("data/nmf.cluster.weight.RData")
### the cluseter weights are determined based on the phylogenetic tree (Ensembl database: http://useast.ensembl.org/index.html)

### function to identify phylogenetic distance of a pair of genes
sr.curr  =sr.all

sr.gene.all = data.table(rescuer = pancancer$genes[sr.curr[,1]], vulnerable = pancancer$genes[sr.curr[,2]])
phylogenetic.screen = phylo.profile.screen(sr.gene.all )
phylogenetic.distance.score = phylo.profile(sr.gene.all )


########## Select pairs that pass Step 4: Phylogenetic screening
sr.phylogenetic = sr.curr[phylogenetic.screen, ]

##########################################################################################
########## Final SR pairs that pass all 4 screenings #####################################
##########################################################################################

####### aggregate stats significance/score of each screen#######

sr.gene.all.screen.stats = sr.gene.all
sr.gene.all.screen.stats$iorio.cohen.score = iorio.cohen.score
sr.gene.all.screen.stats$shRNA.aggreage.fdr = shRNA.significance.fdr
sr.gene.all.screen.stats$sof.aggreage.fdr = sof.aggreage.fdr
sr.gene.all.screen.stats$clinical.aggreage.fdr = clinical.signficance.fdr
sr.gene.all.screen.stats$phylogenetic.distance.score= phylogenetic.distance.score


### final positive SR DU interactions 
positive.sr.du.interactions  = sr.gene.all.screen.stats[which(invitro.screen & clinical.screen & phylogenetic.screen & sof.screen)]
## the data table provide aggregate significance or score of each screenings step 












