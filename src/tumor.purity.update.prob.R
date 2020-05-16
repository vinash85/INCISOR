###update prob with tumor purity 

load("data/TCGA.RData")
library(readxl)
tp =read_excel("~/shortcuts/srescues/ncomms9971-s2_tumor_purity.xlsx",skip=3)
tp = data.table(tp)
setnames(tp,1,"sampleid")
tp[,sample:=gsub(sampleid, pattern = "-01A", replacement = "")]

tcga.tp = tp[match(prob$sample, tp$sample), c("ABSOLUTE", "LUMP", "IHC", "CPE"),with=F]
format.col = function(tt){
	tt=as.numeric(tt)
	ifelse(is.nan(tt), NA,tt)
}
prob$tumor.purity  <- tcga.tp[, lapply(.SD, format.col)]
save(file="data/TCGA.RData", prob)


prob$tumor.purity <- prob$tumor.purity[,lapply(.SD, function(tt) ifelse(is.na(tt), median(tt,na.rm=T),tt))]

pancancer = prob
rm(prob)
### use the SRs ideal number of SR is around 1800 interactions

sr.gene = fread("/cbcb/project2-scratch/jooslee/srescues/pancancer/cox.du/sr.du/sr.final.combined.gene.gii.V2.txt")
sr.all =  cbind(
	match(sr.gene$rescuer, pancancer$genes),
	match(sr.gene$vulnerable, pancancer$genes)
	)




