library(survival)


writeSVMFeature <- function(outfile, pos=NULL, neg=NULL)

{
  outCon = file(outfile, 'w')
  dimLen = dim(pos)[2]
  if(!is.null(pos) && (dimLen > 0) ){
    for(i in seq(1:dim(pos)[1])){
      pos.text = paste("+1",  paste(1:dimLen, pos[i,], sep=":", collapse=" "))
      writeLines(pos.text, con = outCon)
    }
  
  }      
  dimLen = dim(neg)[2]
  if(!is.null(neg) && (dimLen > 0)){
    for(i in seq(1:dim(neg)[1])){
      neg.text = paste("-1",  paste(1:dimLen, neg[i,], sep=":", collapse=" "))
      writeLines(neg.text, con = outCon)
    }        
  }   
  close(outCon)
}


cor.mytest = function(x,y, method1){
	out= rep(NA,2)
	if(sum(!( is.na(x)| is.na(y))) > 2 ){
		xx = cor.test(x,y,method=method1)
		out = c(xx$estimate, xx$p.value)
	}
		out 
}

shrna.screen.achilles = function(sr){
	load("/cbcb/project2-scratch/jooslee/srescues/shrna/prob_essentiality_marcotte.RData")
	ess=prob;rm(prob)
	ess$ess = ess$mat
	sr.ess = cbind( match(pancancer$genes[sr[,1]],ess$genes),
		match(pancancer$genes[sr[,2]],ess$genes))
	require(doMC);require(foreach);registerDoMC(cores = 64)
	outbed = foreach(x = 1:nrow(sr.ess), .inorder=T) %dopar% pair.ess.wil(sr.ess[x,],ess=ess) 
	shrna.wil.mat = do.call(rbind, outbed)
	shrna.wil.mat
}

shrna.screen.bc = function(sr){
	ess.bc = local({load("/cbcb/project2-scratch/jooslee/prob.new/sayad.RData"); environment()})
	ess.bc=ess.bc$dat
	ess.bc$ess = ess.bc$mat
	sr.temp.ess.bc = cbind( match(pancancer$genes[sr[,1]],ess.bc$genes),
		match(pancancer$genes[sr[,2]],ess.bc$genes))
	require(doMC)
	require(foreach)
	registerDoMC(cores = 64)
	outbed = foreach(x = 1:nrow(sr.temp.ess.bc), .inorder=T) %dopar% pair.ess.wil(sr.temp.ess.bc[x,],ess=ess.bc)
	shrna.wil.bc.mat = do.call(rbind, outbed)
	shrna.wil.bc.mat
}


pair.ess.screen = function(pair)
{
		out = rep(NA, 8)
	if(sum(is.na(pair)) == 0){
		a1.mRNA = cor.mytest( ess$mRNA[pair[1],],ess$ess[pair[2],],method="spearman")
		a2.mRNA = cor.mytest( ess$mRNA[pair[2],],ess$ess[pair[1],],method="spearman")
		a1.scna = cor.mytest( ess$scna[pair[1],],ess$ess[pair[2],],method="spearman")
		a2.scna = cor.mytest( ess$scna[pair[2],],ess$ess[pair[1],],method="spearman")
		out =c(a1.mRNA, a2.mRNA, a1.scna, a2.scna)
	}
	out
}
wilcox.test.na = function(a, b , alternative1) {
ifelse(sum(!is.na(a)) > 5 & sum(!is.na(b)) > 5,
	 wilcox.test(a,b,alternative = alternative1)$p.value,
	NA)
}
wilcox.test.na = function(x,y, alternative1, paired=FALSE) {
    tryCatch(
        wilcox.test(x, y, alternative=alternative1, paired=paired)$p.value,
        error = function(e) NA
        )
}

cor.test.na = function(x,y, ...) {
    tryCatch(
        cor.test(x, y, ...),
        error = function(e) list(estimate=NA, p.value=NA)
        )
}
pair.ess.wil = function(pair, ess=ess)
{
		out = rep(NA, 8)
		if(sum(is.na(pair)) == 0){
			aa = ess$mRNA[pair[1],]
			out[1]= wilcox.test.na(ess$ess[pair[2],aa <= median(aa,na.rm=T)],ess$ess[pair[2],aa >median(aa,na.rm=T)],alternative1="less" )
			aa = ess$mRNA[pair[2],]
			out[2]= wilcox.test.na(ess$ess[pair[1],aa <= median(aa,na.rm=T)],ess$ess[pair[1],aa >median(aa,na.rm=T)],alternative1="less" )
			aa = ess$scna[pair[1],]
			out[3]= wilcox.test.na(ess$ess[pair[2],aa <= median(aa,na.rm=T)],ess$ess[pair[2],aa >median(aa,na.rm=T)],alternative1="less" )
			aa = ess$scna[pair[2],]
			out[4]= wilcox.test.na(ess$ess[pair[1],aa <= median(aa,na.rm=T)],ess$ess[pair[1],aa >median(aa,na.rm=T)],alternative1="less" )

			aa = ess$mRNA[pair[1],]
			out[5]= wilcox.test.na(ess$ess[pair[2],aa <= median(aa,na.rm=T)],ess$ess[pair[2],aa >median(aa,na.rm=T)],alternative1="greater" )
			aa = ess$mRNA[pair[2],]
			out[6]= wilcox.test.na(ess$ess[pair[1],aa <= median(aa,na.rm=T)],ess$ess[pair[1],aa >median(aa,na.rm=T)],alternative1="greater" )
			aa = ess$scna[pair[1],]
			out[7]= wilcox.test.na(ess$ess[pair[2],aa <= median(aa,na.rm=T)],ess$ess[pair[2],aa >median(aa,na.rm=T)],alternative1="greater" )
			aa = ess$scna[pair[2],]
			out[8]= wilcox.test.na(ess$ess[pair[1],aa <= median(aa,na.rm=T)],ess$ess[pair[1],aa >median(aa,na.rm=T)],alternative1="greater" )
	}
	out
}

pair.ccle.wil = function(target, ccle=ccle, probs = .5, alternative="less"){
		out = matrix(NA, 2*nrow(ccle$mat))
		require(parallel)

			aa = ccle$mRNA[target,]
			thr = quantile(aa,probs=probs, na.rm=T)
			out[1:nrow(ccle$mat)]= sapply(seq(nrow(ccle$mat)), function(tt) {
				wilcox.test.na(ccle$mat[tt,aa <= thr ], ccle$mat[tt,aa >thr],alternative1= alternative ) ## cell  line lower expression of genes have lower IC50 and hence more effective 
			})
			aa = ccle$scna[target,]
			thr = quantile(aa,probs=probs, na.rm=T)
			out[(1+nrow(ccle$mat)) : (2*nrow(ccle$mat))] = sapply(seq(nrow(ccle$mat)), function(tt) {
				wilcox.test.na(ccle$mat[tt,aa <= thr ], ccle$mat[tt,aa >thr],alternative1=alternative )
				}) 

	out
}



mol.pattern = function(inp){
	zero.inx = inp< 1E-300
	one.inx = inp ==1 
	inp[inp < 0.5] = -inp[inp< 0.5]
	inp[inp > 0.5] = 1 - inp[inp > 0.5]
	inp[zero.inx] = -1E-300
	inp[one.inx] = 1E-300

	inp = -inp
	inp[c(3,6)] = -inp[c(3,6)]
	inp = ifelse(inp<0, 1, inp) 
	t1 = max(min(inp[1:2]),inp[3]) 
	t2 = max(min(inp[4:5]),inp[6])
	out = ifelse(t1< 0.1 & t2 < 0.1 , t1*t2, 1 )
	# t2 = ifelse(t2< 1, t2, 1 )
	out 
}


cox.sl = function(pairs, scnaq, scna, survival){
	out = matrix(0, ncol=2, nrow=nrow(pairs))
	for (ii in 1:nrow(pairs)) {
		pair = pairs[ii,]
		g1 = scna[pair[1],]  
		g2 = scna[pair[2],]
		f1 = scnaq[pair[1],]  
		f2 = scnaq[pair[2],]  
		cov = ifelse(f1 == 0 & f2==0, 1, 0 )
		dt1 = data.table(survival, g1 , g2, cov=cov)
		setnames(dt1,1:2,c("time", "status"))
		cox.out = coxph(Surv(time,status) ~., data=dt1)
		aa  = summary(cox.out)
		bb = aa$coefficients["cov",]
		pval = ifelse(bb[1] < 0 , bb[5],1)
		out[ii,] = c(bb[2], pval)
	}
	out
} 
cox.sdl = function(pairs, scnaq, scna, survival){
	out = matrix(0, ncol=2, nrow=nrow(pairs))
	for (ii in 1:nrow(pairs)) {
		pair = pairs[ii,]
		g1 = scna[pair[1],]  
		g2 = scna[pair[2],]
		f1 = scnaq[pair[1],]  
		f2 = scnaq[pair[2],]  
		cov = ifelse(f1 == 0 & f2==0, 1, 0 )
		dt1 = data.table(survival, g1 , g2, cov=cov)
		setnames(dt1,1:2,c("time", "status"))
		cox.out = coxph(Surv(time,status) ~., data=dt1)
		aa  = summary(cox.out)
		bb = aa$coefficients["cov",]
		pval = ifelse(bb[1] < 0 , bb[5],1)
		out[ii,] = c(bb[2], pval)
	}
	out
} 
cox.pair.sl = function(pair,use.mRNA=F)
{
	if(use.mRNA){
		g1 = mRNA.norm[pair[1],]  
		g2 = mRNA.norm[pair[2],]
		f1 = mRNAq2[pair[1],]  
		f2 = mRNAq2[pair[2],]  
	}else{
		g1 = scna.norm[pair[1],]  
		g2 = scna.norm[pair[2],]  
		f1 = scnaq2[pair[1],]  
		f2 = scnaq2[pair[2],]  
	}
	dt1 = cbind(surv.type, cbind(g1 , g2))
	cntrl.out = coxph(Surv(time,status) ~., data=dt1)
	ll1 = cntrl.out$loglik[2]
	cov = ifelse(f1 == 0 & f2==0, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	dt2 = surv.type
	dt2$cov = qnorm.array(cov)
	cox.out.nocntrl = coxph(Surv(time,status) ~., data=dt2)
	cc  = summary(cox.out.nocntrl)
	c(aa$coefficients["cov",], ll2 - ll1, cc$coefficients["cov",])
}

cox.pair.sdl = function(pair,use.mRNA=F)
{
	if(use.mRNA){
		g1 = mRNA.norm[pair[1],]  
		g2 = mRNA.norm[pair[2],]
		f1 = pancancer$mRNAq2[pair[1],]  
		f2 = pancancer$mRNAq2[pair[2],]  
	}else{
		g1 = scna.norm[pair[1],]  
		g2 = scna.norm[pair[2],]  
		f1 = pancancer$scnaq2[pair[1],]  
		f2 = pancancer$scnaq2[pair[2],]  
	}
	cov = ifelse(f1 == 0 & f2 == 2,1,0 )
	dt1 = cbind(surv.type, cbind(g1 , g2))
	cntrl.out = coxph(Surv(time,status) ~., data=dt1)
	ll1 = cntrl.out$loglik[2]
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	dt1 = surv.type
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	bb  = summary(cox.out)
	
	c(aa$coefficients["cov",], ll2 - ll1, bb$coefficients["cov",])
} 


cox.pair.du = function(pair,use.mRNA=F)
{
	if(use.mRNA){
		g1 = mRNA.norm[pair[1],]  
		g2 = mRNA.norm[pair[2],]
		f1 = mRNAq2[pair[1],]  
		f2 = mRNAq2[pair[2],]  
	}else{
		g1 = scna.norm[pair[1],]  
		g2 = scna.norm[pair[2],]  
		f1 = scnaq2[pair[1],]  
		f2 = scnaq2[pair[2],]  
	}
	# is.nan.inx = which((!(is.na(g1))) & (!(is.na(g1))) )
	dt1 = cbind(surv.type, cbind(g1 , g2))
	cntrl.out = coxph(Surv(time,status) ~., data=dt1)
	ll1 = cntrl.out$loglik[2]
	cov = ifelse(f2 == 0 & f1<=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 0 & f1==2, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	uu = c(aa$coefficients["cov",], ll2 - ll1, bb$coefficients["cov",], ll3-ll1)
	dt1 = surv.type
	cov = ifelse(f2 == 0 & f1<=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 0 & f1==2, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	c(uu, aa$coefficients["cov",],  bb$coefficients["cov",])
	
} 

cox.pair.dd = function(pair,use.mRNA=F)
{
	if(use.mRNA){
		g1 = mRNA.norm[pair[1],]  
		g2 = mRNA.norm[pair[2],]
		f1 = mRNAq2[pair[1],]  
		f2 = mRNAq2[pair[2],]  
	}else{
		g1 = scna.norm[pair[1],]  
		g2 = scna.norm[pair[2],]  
		f1 = scnaq2[pair[1],]  
		f2 = scnaq2[pair[2],]  
	}
	dt1 = cbind(surv.type, cbind(g1 , g2))
	cntrl.out = coxph(Surv(time,status) ~., data=dt1)
	ll1 = cntrl.out$loglik[2]
	cov = ifelse(f2 == 0 & f1>=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 0 & f1==0, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	uu = c(aa$coefficients["cov",], ll2 - ll1, bb$coefficients["cov",], ll3-ll1)
	dt1 = surv.type
	cov = ifelse(f2 == 0 & f1>=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 0 & f1==0, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	c(uu, aa$coefficients["cov",],  bb$coefficients["cov",])
	
}

cox.pair.ud = function(pair,use.mRNA=F)
{
	if(use.mRNA){
		g1 = mRNA.norm[pair[1],]  
		g2 = mRNA.norm[pair[2],]
		f1 = mRNAq2[pair[1],]  
		f2 = mRNAq2[pair[2],]  
	}else{
		g1 = scna.norm[pair[1],]  
		g2 = scna.norm[pair[2],]  
		f1 = scnaq2[pair[1],]  
		f2 = scnaq2[pair[2],]  
	}
	dt1 = cbind(surv.type, cbind(g1 , g2))
	cntrl.out = coxph(Surv(time,status) ~., data=dt1)
	ll1 = cntrl.out$loglik[2]
	cov = ifelse(f2 == 2 & f1>=1, 1, 0 )

	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 2 & f1==0, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	uu = c(aa$coefficients["cov",], ll2 - ll1, bb$coefficients["cov",], ll3-ll1)
	dt1 = surv.type
	cov = ifelse(f2 == 2 & f1>=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 2 & f1==0, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	c(uu, aa$coefficients["cov",],  bb$coefficients["cov",])
	
}




cox.pair.uu = function(pair,use.mRNA=F)
{
	if(use.mRNA){
		g1 = mRNA.norm[pair[1],]  
		g2 = mRNA.norm[pair[2],]
		f1 = mRNAq2[pair[1],]  
		f2 = mRNAq2[pair[2],]  
	}else{
		g1 = scna.norm[pair[1],]  
		g2 = scna.norm[pair[2],]  
		f1 = scnaq2[pair[1],]  
		f2 = scnaq2[pair[2],]  
	}
	# is.nan.inx = which((!(is.na(g1))) & (!(is.na(g1))) )
	dt1 = cbind(surv.type, cbind(g1 , g2))
	cntrl.out = coxph(Surv(time,status) ~., data=dt1)
	ll1 = cntrl.out$loglik[2]
	cov = ifelse(f2 == 2 & f1<=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 2 & f1==2, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	uu = c(aa$coefficients["cov",], ll2 - ll1, bb$coefficients["cov",], ll3-ll1)
	dt1 = surv.type
	cov = ifelse(f2 == 2 & f1<=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 2 & f1==2, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	c(uu, aa$coefficients["cov",],  bb$coefficients["cov",])
	
} 

cox.pair.sl.strata = function(pair,use.mRNA=F)
{
	if(use.mRNA){
		g1 = mRNA.norm[pair[1],]  
		g2 = mRNA.norm[pair[2],]
		f1 = mRNAq2[pair[1],]  
		f2 = mRNAq2[pair[2],]  
	}else{
		g1 = scna.norm[pair[1],]  
		g2 = scna.norm[pair[2],]  
		f1 = scnaq2[pair[1],]  
		f2 = scnaq2[pair[2],]  
	}
	dt1 = cbind(surv.strata, cbind(g1 , g2))
	cntrl.out = coxph(Surv(time,status) ~ g1 + g2 + strata(type), data=dt1)
	ll1 = cntrl.out$loglik[2]
	cov = ifelse(f1 == 0 & f2==0, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~ g1 + g2 + cov + strata(type), data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	dt2 = surv.strata
	dt2$cov = qnorm.array(cov)
	cox.out.nocntrl = coxph(Surv(time,status) ~., data=dt2)
	cc  = summary(cox.out.nocntrl)
	c(aa$coefficients["cov",], ll2 - ll1, cc$coefficients["cov",] )
} 

cox.pair.du.strata = function(pair,use.mRNA=F)
{
	if(use.mRNA){
		g1 = mRNA.norm[pair[1],]  
		g2 = mRNA.norm[pair[2],]
		f1 = mRNAq2[pair[1],]  
		f2 = mRNAq2[pair[2],]  
	}else{
		g1 = scna.norm[pair[1],]  
		g2 = scna.norm[pair[2],]  
		f1 = scnaq2[pair[1],]  
		f2 = scnaq2[pair[2],]  
	}
	dt1 = cbind(surv.strata, cbind(g1 , g2))
	cntrl.out = coxph(Surv(time,status) ~ g1 + g2 + strata(type), data=dt1)
	ll1 = cntrl.out$loglik[2]
	cov = ifelse(f2 == 0 & f1<=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~ g1 + g2 + cov + strata(type), data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 0 & f1==2, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~ g1 + g2 + cov + strata(type) , data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	uu = c(aa$coefficients["cov",], ll2 - ll1, bb$coefficients["cov",], ll3-ll1)
	dt1 = surv.strata
	cov = ifelse(f2 == 0 & f1<=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~  cov + strata(type), data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 0 & f1==2, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~ cov + strata(type), data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	c(uu, aa$coefficients["cov",],  bb$coefficients["cov",])
	
} 






cox.pair.du.strata.controlled = function(pair,prob, use.mRNA=F)
{

	coxph.run = function(formula.list,data1){
		fmla <- as.formula(paste(" Surv(time,status) ~ ", paste(formula.list, collapse= "+")))
		cox.out = coxph(fmla, data=data1)
		ll = cox.out$loglik[2]
		aa  = summary(cox.out)
		out = ll 
		if("cov" %in% rownames(aa$coefficients) ){
			out = c(aa$coefficients["cov", ], ll)
		}
		out
	}
	coxph.robust = function(data1, f1.list, f2.list=NULL){
		tryCatch( 
			coxph.run(c(f1.list, f2.list), data1),
			error = function(e)  coxph.run(f1.list, data1)
			)
	}


	if(use.mRNA){
		g1 = prob$mRNA.norm[pair[1],]  
		g2 = prob$mRNA.norm[pair[2],]
		f1 = prob$mRNAq2[pair[1],]  
		f2 = prob$mRNAq2[pair[2],]  
	}else{
		g1 = prob$scna.norm[pair[1],]  
		g2 = prob$scna.norm[pair[2],]  
		f1 = prob$scnaq2[pair[1],]  
		f2 = prob$scnaq2[pair[2],]  
	}
	surv.strata = prob$surv.strata
	f2.list = c("strata(types)",  "strata(sex)",  "age", "strata(race)") 
	dt1 = cbind(surv.strata, cbind(g1 , g2))
	cntrl.out = coxph.robust(f1.list = c("g1",  "g2",  "strata(type)"), f2.list = f2.list,data1=dt1)
	ll1 = cntrl.out
	cov = ifelse(f2 == 0 & f1<=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph.robust(f1.list = c("g1",  "g2",  "strata(type)", "cov"), f2.list = f2.list,data1=dt1) 
	ll2 = cox.out[6]
	aa  = cox.out[1:5]
	cov = ifelse(f2 == 0 & f1==2, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph.robust(f1.list = c("g1",  "g2",  "strata(type)", "cov"), f2.list = f2.list,data1=dt1) 
	ll3 = cox.out[6]
	bb  = cox.out[1:5]
	uu = c(aa, ll2 - ll1, bb, ll3-ll1)
	
	dt1 = surv.strata
	cov = ifelse(f2 == 0 & f1<=1, 1, 0 )
	dt1$cov = qnorm.array(cov)


	cox.out = coxph.robust(f1.list = c(  "strata(type)", "cov"), f2.list = f2.list,data1=dt1)  
	aa  = cox.out[1:5]
	cov = ifelse(f2 == 0 & f1==2, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph.robust(f1.list = c(  "strata(type)", "cov"), f2.list = f2.list,data1=dt1)  
	bb  = cox.out[1:5]
	c(uu, aa,  bb)
	
} 


cox.pair.dd.strata = function(pair,use.mRNA=F)
{
	if(use.mRNA){
		g1 = mRNA.norm[pair[1],]  
		g2 = mRNA.norm[pair[2],]
		f1 = mRNAq2[pair[1],]  
		f2 = mRNAq2[pair[2],]  
	}else{
		g1 = scna.norm[pair[1],]  
		g2 = scna.norm[pair[2],]  
		f1 = scnaq2[pair[1],]  
		f2 = scnaq2[pair[2],]  
	}
	dt1 = cbind(surv.strata, cbind(g1 , g2))
	cntrl.out = coxph(Surv(time,status) ~ g1 + g2 + strata(type), data=dt1)
	ll1 = cntrl.out$loglik[2]
	cov = ifelse(f2 == 0 & f1>=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~ g1 + g2 + cov + strata(type), data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 0 & f1==0 , 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~ g1 + g2 + cov + strata(type) , data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	uu = c(aa$coefficients["cov",], ll2 - ll1, bb$coefficients["cov",], ll3-ll1)
	dt1 = surv.strata
	cov = ifelse(f2 == 0 & f1>=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~  cov + strata(type), data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 0 & f1==0, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~  cov + strata(type), data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	c(uu, aa$coefficients["cov",],  bb$coefficients["cov",])
	
} 

cox.pair.ud.strata = function(pair,use.mRNA=F)
{
	if(use.mRNA){
		g1 = mRNA.norm[pair[1],]  
		g2 = mRNA.norm[pair[2],]
		f1 = mRNAq2[pair[1],]  
		f2 = mRNAq2[pair[2],]  
	}else{
		g1 = scna.norm[pair[1],]  
		g2 = scna.norm[pair[2],]  
		f1 = scnaq2[pair[1],]  
		f2 = scnaq2[pair[2],]  
	}
	dt1 = cbind(surv.strata, cbind(g1 , g2))
	cntrl.out = coxph(Surv(time,status) ~ g1 + g2 + strata(type), data=dt1)
	ll1 = cntrl.out$loglik[2]
	cov = ifelse(f2 == 2 & f1>=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~ g1 + g2 + cov + strata(type), data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 2 & f1==0 , 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~ g1 + g2 + cov + strata(type) , data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	uu = c(aa$coefficients["cov",], ll2 - ll1, bb$coefficients["cov",], ll3-ll1)
	dt1 = surv.strata
	cov = ifelse(f2 == 2 & f1>=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~  cov + strata(type), data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 2 & f1==0, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~  cov + strata(type), data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	c(uu, aa$coefficients["cov",],  bb$coefficients["cov",])
	
} 

cox.pair.uu.strata = function(pair,use.mRNA=F)
{
	if(use.mRNA){
		g1 = mRNA.norm[pair[1],]  
		g2 = mRNA.norm[pair[2],]
		f1 = mRNAq2[pair[1],]  
		f2 = mRNAq2[pair[2],]  
	}else{
		g1 = scna.norm[pair[1],]  
		g2 = scna.norm[pair[2],]  
		f1 = scnaq2[pair[1],]  
		f2 = scnaq2[pair[2],]  
	}
	dt1 = cbind(surv.strata, cbind(g1 , g2))
	cntrl.out = coxph(Surv(time,status) ~ g1 + g2 + strata(type), data=dt1)
	ll1 = cntrl.out$loglik[2]
	cov = ifelse(f2 == 2 & f1<=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~ g1 + g2 + cov + strata(type), data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 2 & f1==2, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~ g1 + g2 + cov + strata(type) , data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	uu = c(aa$coefficients["cov",], ll2 - ll1, bb$coefficients["cov",], ll3-ll1)
	dt1 = surv.strata
	cov = ifelse(f2 == 2 & f1<=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~  cov + strata(type), data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 2 & f1==2, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~ cov + strata(type), data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	c(uu, aa$coefficients["cov",],  bb$coefficients["cov",])
	
} 


cox.pair.sl.type = function(pair,use.mRNA=F){
	if(use.mRNA){
		g1 = mRNA.norm.sub[pair[1],]  
		g2 = mRNA.norm.sub[pair[2],]
		f1 = mRNAq2.sub[pair[1],]  
		f2 = mRNAq2.sub[pair[2],]  
	}else{
		g1 = scna.norm.sub[pair[1],]  
		g2 = scna.norm.sub[pair[2],]  
		f1 = scnaq2.sub[pair[1],]  
		f2 = scnaq2.sub[pair[2],]  
	}

	dt1 = cbind(surv.type.dt, cbind(g1 , g2))
	cntrl.out = coxph(Surv(time,status) ~., data=dt1)
	ll1 = cntrl.out$loglik[2]
	cov = ifelse(f2 == 0 & f1==0, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	# uu = c(aa$coefficients["cov",], ll2 - ll1, bb$coefficients["cov",], ll3-ll1)
	dt1 = surv.type.dt
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	bb  = summary(cox.out)
	c( aa$coefficients["cov",], ll2 - ll1,   bb$coefficients["cov",])
} 



cox.pair.du.type = function(pair,use.mRNA=F){
	if(use.mRNA){
		g1 = mRNA.norm.sub[pair[1],]  
		g2 = mRNA.norm.sub[pair[2],]
		f1 = mRNAq2.sub[pair[1],]  
		f2 = mRNAq2.sub[pair[2],]  
	}else{
		g1 = scna.norm.sub[pair[1],]  
		g2 = scna.norm.sub[pair[2],]  
		f1 = scnaq2.sub[pair[1],]  
		f2 = scnaq2.sub[pair[2],]  
	}

	dt1 = cbind(surv.type.dt, cbind(g1 , g2))
	cntrl.out = coxph(Surv(time,status) ~., data=dt1)
	ll1 = cntrl.out$loglik[2]
	cov = ifelse(f2 == 0 & f1<=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 0 & f1==2, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	uu = c(aa$coefficients["cov",], ll2 - ll1, bb$coefficients["cov",], ll3-ll1)
	dt1 = surv.type.dt
	cov = ifelse(f2 == 0 & f1<=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 0 & f1==2, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	c(uu, aa$coefficients["cov",],  bb$coefficients["cov",])
} 

cox.pair.dd.type = function(pair,use.mRNA=F){
	if(use.mRNA){
		g1 = mRNA.norm.sub[pair[1],]  
		g2 = mRNA.norm.sub[pair[2],]
		f1 = mRNAq2.sub[pair[1],]  
		f2 = mRNAq2.sub[pair[2],]  
	}else{
		g1 = scna.norm.sub[pair[1],]  
		g2 = scna.norm.sub[pair[2],]  
		f1 = scnaq2.sub[pair[1],]  
		f2 = scnaq2.sub[pair[2],]  
	}
	dt1 = cbind(surv.type.dt, cbind(g1 , g2))
	cntrl.out = coxph(Surv(time,status) ~., data=dt1)
	ll1 = cntrl.out$loglik[2]
	cov = ifelse(f2 == 0 & f1>=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 0 & f1==0, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	uu = c(aa$coefficients["cov",], ll2 - ll1, bb$coefficients["cov",], ll3-ll1)	
	
	dt1 = surv.type.dt
	cov = ifelse(f2 == 0 & f1>=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	cov = ifelse(f2 == 0 & f1==0, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph(Surv(time,status) ~., data=dt1)
	ll3 = cox.out$loglik[2]
	bb  = summary(cox.out)
	c(uu, aa$coefficients["cov",],  bb$coefficients["cov",])
} 



 

co.exp = function(pair){
	g1 = pancancer$mRNA[pair[1],]  
	g2 = pancancer$mRNA[pair[2],]  
	aa = cor.test(g1, g2, method = "spearman")
	g1 = pancancer$scna[pair[1],]  
	g2 = pancancer$scna[pair[2],]  
	bb = cor.test(g1, g2, method = "spearman")
	c(aa$estimate, aa$p.value, bb$estimate, bb$p.value)
}


create.cd.fig = function(sr.curr, vDelcells, vNormcells, scnaqtest, survType, weights = NULL, type="sr", labels=list(Del="not-rescued", Norm="rescued"), ppt=F){
	if(is.null(weights)) weights= rep(1, nrow(sr.curr))
	survivaltest = survType[,1:2]
	vNorm = vLet = vDel=numeric(ncol(scnaqtest))
	for (idx in 1:nrow(sr.curr)){
		gene1 = sr.curr[idx, 1]
		gene2 = sr.curr[idx, 2]
		scna2 = scnaqtest[gene2,]
		scna1 = scnaqtest[gene1, ]
		inx = 3 * scna2 + scna1 + 1;
		vDel= vDel+ (inx %in% vDelcells) * weights[idx]
		vNorm= vNorm + (inx %in%  vNormcells);
	}
	vInx = which(vNorm > 0)
	vDel.s = vDel[vInx]
	vNormT=vDel.s<=quantile(vDel.s,0.1)
	vDelT=vDel.s>=quantile(vDel.s,0.9)
	# browser()
	times1=survivaltest[vInx[vDelT],,drop=F];
	times2=survivaltest[vInx[vNormT],,drop=F];
	outf3= logRank(times1,times2);
	###
	require(data.table)
	require(survival)
	require(ggplot2)
	dt = data.table( rbind(times1, times2)  )
	if(type=="sl"){
		dt$label = "SL-"
		dt$label[1:nrow(times1)] = "SL+"
	}else{
		dt$label = labels$Del
		dt$label[1:nrow(times1)] = labels$Norm
	}
	setnames(dt, 1:2, c("time", "status"))
	source("/cbcb/project2-scratch/jooslee/srescues/shrna/ggsurv.R")
	sr.surv <- survfit(Surv(time,status==0)~label, data = dt)
	# browser() 
	surv.dt1 = data.table(survType)
	surv.dt1$cov = vDel
	setnames(surv.dt1, 1:2, c("time", "status"))
	# cox.out = coxph(Surv(time,status==0) ~., data=surv.dt1)
	cox.out = coxph(Surv(time,status==0) ~., data=surv.dt1)
	bb  = summary(cox.out)
	 cc = bb$coefficients["cov",]

	dt$marker <- c(rep(0,nrow(times1)), rep(1,nrow(times2)))
	dt$status1 = ifelse(dt$status==0,1,0)
	aa <- ggsurv(sr.surv)
	text1 = paste0("P=", formatC(outf3[1], format="E", digits=2), "\n")
	text2 = paste0("AUC=", formatC(outf3[8] - outf3[7],digits=2))
	text = paste0(text1, " delta ", text2 )

	coxpx = formatC(cc[5], format="E", digits=2)
	betax = formatC(cc[2], format="E", digits=2)
	px = formatC(outf3[1], format="E", digits=2)
	aucx = formatC(outf3[7] - outf3[8],digits=2)
	aa <- ggsurv(sr.surv)
	text = paste0(
                "\u0394 AUC", "=", aucx, "\n", "P=", px, "\n",
                "Hazard ratio", "=",  betax, "\n", "Cox-P=", coxpx 
                
                ) 
	if(ppt){

		aa <-aa +
		annotate( "text", x = 1000, y = .4, label = text) 
	}else{

		aa <-aa +
		annotate( "text", x = 2000, y = .4,  parse=F, label = as.character(text)) +
		theme(axis.title=element_text(size=20), legend.text=element_text(size=20)) +
		theme(strip.background=element_blank(), strip.text.y = element_text()) +
		theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank())  
	}
	
aa
}

qnorm.array <- function(mat)
{
	mat.back = mat 
	mat = mat.back[!is.na(mat.back)]
    mat = rank(mat,  rank, ties.method = "average");
    mat = qnorm(mat / (length(mat)+1));
    mat.back[!is.na(mat.back)] = mat 
    mat.back
}


cox.single.gene.strata = function(inx, g1, surv.strata)
{
	dt1 = cbind(surv.strata, g1)
	cntrl.out = coxph(Surv(time,status) ~ g1 +  strata(types), data=dt1)
	ll1 = cntrl.out$loglik[2]
	aa  = summary(cntrl.out)
	c(aa$coefficients["g1",], ll1) 
} 



phylo.profile = function(sr.gene.all){
	load("/cbcb/project2-scratch/jooslee/srescues/pancancer/cox.du/yuval.phylogenetic.profile.RData")
	load("/cbcb/project2-scratch/jooslee/srescues/pancancer/cox.du/feature.weight.RData")
	sr.gene1 = sr.gene.all
	sr.phylo =  cbind(match(sr.gene1$rescuer, phylo$genes), match(sr.gene1$vulnerable, phylo$genes))
	featureMat = (phylo[sr.phylo[,1],-(1:3)] - phylo[sr.phylo[,2],-(1:3)])^2
	featureMat %*% t(feature.weight)
}


rank.norm= function(ss) {
	out = rank(ss, na.last="keep")/max(1, sum(!is.na(ss)))
	out[is.na(out)] = 1
	out
}
