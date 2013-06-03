#######################################################
# Function: robust.fdr
# Purpose: Implement robust method of FDR estimation (Pounds and Cheng 2006, Bioinformatics)
# Arguments: p - vector of p-values from the analysis
#            sides - indicate whether p-values are 1-sided (set sides=1) or 2-sided (set sides=2), default=1
#            p2 - for one-sided testing, p-values from testing the "other alternative", default=1-p
#            discrete - indicates (T/F) whether p-values are discrete
#            use8 - indicates whether the constant 8 should be used if p-values are discrete, see Pounds and Cheng (2006) for more details
# Returns: a list with the following components
#            p - the vector of p-values provided by the user
#            fdr - the vector of smoothed FDR estimates
#            q - the vector of q-values based on the smoothed FDR estimates
#            cdf - the vector with p-value empirical distribution function at corresponding entry of p
#            loc.fdr - the local (unsmoothed) FDR estimates
#            fp - the estimated number of false positives at p-value cutoff in p
#            fn - the estimated number of false negatives at p-value cutoff in p
#            te - the total of fp and fn
#            pi - the null proportion estimate
#            ord - a vector of indices to order the vectors above by ascending p-value
#########################################################

robust.fdr<-function(p,sides=1,p2=1-p,discrete=F,use8=T)

{
  require(MASS)
  lts.rank.regr<-function(x,y){
    rx<-rank(x)
    ry<-rank(y)
    rfit<-ltsreg(rx,ry)
    ryhat<-rfit$fitted.values
    yhat<-approx(ry,y,xout=ryhat,rule=2)$y
    return(list(x=x,y=y,yhat=yhat))
  }
	m<-length(p)
	ord<-order(p)
	pcumtab<-cumsum(table(p))/m
	F05<-mean(p<=0.5)	
	edf<-approx(as.numeric(names(pcumtab)),pcumtab,xout=p,rule=2)$y
	if (sides==2){
		pi<-min(1,2*mean(p))
		loc.fdr<-pi*p/edf
	}
	else{
		p.new<-2*(p*(p<=p2)+p2*(p>p2))
		pi<-min(1,2*mean(p.new))
		if (discrete) {
			if (use8) pi<-min(1,8*mean(p.new))
			else{
				lam<-max(p[p<=0.5])
				k<-1/(lam^2+(0.5-lam)^2)
				pi<-min(k*mean(p.new),1)
			
			}
		}
		loc.fdr<-pi*p/edf
		loc.fdr[p>0.5]<-(0.5*pi+edf[p>0.5]-F05)/edf[p>0.5]
	}
	up<-unique(p)
	ufdr<-approx(p[ord],loc.fdr[ord],xout=up)$y
	res<-lts.rank.regr(up,ufdr)
	fdr<-approx(up,res$yhat,xout=p)$y
	fdr[fdr>1]<-1
	fp<-m*edf*fdr
	fn<-m*(F05-edf-pi*(0.5-p))*(p<0.5)*(sides==1)+(sides==2)*m*(1-edf-pi*(1-p))
	fn[fn<0]<-0
	te<-fp+fn
	q<-fdr
	q[rev(ord)]<-cummin(fdr[rev(ord)])
	
	return(list(p=p,fdr=fdr,q=q,cdf=edf,loc.fdr=loc.fdr,fp=fp,fn=fn,te=te,pi=pi,ord=ord))
}