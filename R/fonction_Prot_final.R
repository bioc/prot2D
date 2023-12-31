require(fdrtool)
require(st)
require(Biobase)
require(limma)
require(samr)
require(qvalue)
require(Mulcom)

############ Normalize the data and visualize normalization ############

Norm.qt <- 
  function(data, n1, n2, plot=T){
    if(class(data)!="data.frame" & class(data)!="matrix"){
      warning(gettextf("'%s' is neither a dataframe nor a matrix",
                       as.character(match.call()$data),domain=NA))
    }
    if(ncol(data)!=n1+n2){
      warning(gettextf("The number of columns in '%s' doesn't match the number of replicates ", 
                       as.character(match.call()$data),domain=NA))
    }
    mat <- as.matrix(data)
    mat.qt <- normalizeQuantiles(mat)
    
    if(plot){
      par(mfrow=c(2,1))
      RIplot(data,n1=n1,n2=n2, main="Raw data")
      RIplot(mat.qt,n1=n1,n2=n2, main="Normalized data")
      par(mfrow=c(1,1))
    }
    norm <- log2(mat.qt)
    colnames(norm)<- colnames(data)
    rownames(norm) <- rownames(data)
    return(norm)
  }

Norm.vsn <- 
  function(data, n1, n2, plot=T){
    if(class(data)!="data.frame" & class(data)!="matrix"){
      warning(gettextf("'%s' is neither a dataframe nor a matrix",
                       as.character(match.call()$data),domain=NA))
    }
    if(ncol(data)!=n1+n2){
      warning(gettextf("The number of columns in '%s' doesn't match the number of replicates ", 
                       as.character(match.call()$data),domain=NA))
    }
    mat <- as.matrix(data)
    mat.vsn <- normalizeVSN(mat, verbose=F)
    
    if(plot){
      par(mfrow=c(2,1))
      RIplot(data,n1=n1,n2=n2, main="Raw data")
      RIplot(2^mat.vsn,n1=n1,n2=n2,main="Normalized data")
      par(mfrow=c(1,1))
    }
    colnames(mat.vsn)<- colnames(data)
    rownames(mat.vsn) <- rownames(data)
    return(mat.vsn)
  }  

#### Function for R-I plot ####

RIplot<-
  function(data, n1, n2,...){
    if(class(data)!="data.frame" & class(data)!="matrix"){
      warning(gettextf("'%s' is neither a dataframe nor a matrix",
                       as.character(match.call()$data),domain=NA))
    }
    if(ncol(data)!=n1+n2){
      warning(gettextf("The number of columns in '%s' doesn't match the number of replicates ", 
                       as.character(match.call()$data),domain=NA))
    }
    c1.r <- data[,1:n1]
    c2.r <- data [,(n1+1):(n1+n2)]
    m.c1.r <- apply(c1.r, 1, mean)
    m.c2.r <- apply (c2.r,1,mean)
    ratio.r <- log2(m.c2.r/m.c1.r)
    int.r <- log10(m.c1.r*m.c2.r)
    subset(data, ratio.r>1)-> r.UP
    subset(data, ratio.r< -1)-> r.DO
    plot(ratio.r~int.r, pch=20, col=grey(0.5), 
         xlab="Intensity", ylab="Ratio",...)
    abline(h=c(-1,0,1))
    points(ratio.r[ratio.r>1]~int.r[ratio.r>1], pch=20, col="#EE6A50")
    points(ratio.r[ratio.r< -1]~int.r[ratio.r< -1], pch=20, col="#8EE5EE")
    text (x=max(int.r)-1, y=1.5, labels=dim(r.UP)[1])
    text (x=max(int.r)-1, y= -1.5, labels=dim(r.DO)[1])
    res <- data.frame(ratio=ratio.r, intensity =int.r)
    invisible(list(RI=res, diff=c(dim(r.UP)[1],dim(r.DO)[1])))
  }

##### Fonction for gouping data in an Expression Set ###
### data = matrix spots intensities with gel as colums (with replicates of one condition
# following) and spots as raws
# normally this returns by the "Norm.loess" function
#data are compute in the "assayData" slot (and log2 normalized)
#### factor= factor should be in a form of a dataframe with 1 column giving the 2 levels
# for the gels and the names of gels as rownames (in a same order as data!)
#Compute in the slot "featureData" the ratio (log2)
ES.prot <- 
  function(data,n1,n2, f){
    if(class(data)!="matrix"){
      warning(gettextf("'%s' is not a matrix",
                       as.character(match.call()$data),domain=NA))
    }
    if(class(f)!="data.frame"){
      warning(gettextf("'%s' is not a dataframe",
                       as.character(match.call()$f),domain=NA))
    }
    if(ncol(data)!=n1+n2){
      warning(gettextf("The number of columns in '%s' doesn't match the number of replicates ",
                       as.character(match.call()$data),domain=NA))  
    }
    if(is.null(colnames(data))){
      warning(gettextf("Columns of '%s' don't have names",
                       as.character(match.call()$data),domain=NA))  
    }
    if(any(colnames(data)!=rownames(f))){
      warning(gettextf("The names of columns in '%1$s' doesn't match the names of rows in '%2$s' ",
                       as.character(match.call()$data),as.character(match.call()$f),domain=NA))          
      
    }
    c1 <- data[,1:n1]
    c2 <- data[,(n1+1):(n1+n2)]
    m.c1 <- apply(c1, 1, mean)
    m.c2 <- apply (c2,1,mean)
    ratio <-m.c2- m.c1
    r <- new("AnnotatedDataFrame", data=as.data.frame(ratio))
    f <- new("AnnotatedDataFrame", data=f)
    ES <-ExpressionSet(assayData=data, phenoData=f, featureData=r)
    return(ES)
  }


######## Functions for finding differentially expressed spots #######

#### Student's T-test #######
ttest.Prot <-
  function(data, plot = T, fdr.thr = 0.1, Fold2 = F, method.fdr = "BH", var.equal = F){
    if(class(data)[1]!="ExpressionSet"){
      warning(gettextf("'%s' is not an ExpressionSet",
                       as.character(match.call()$data),domain=NA))  
    }
    if(method.fdr != "BH" & method.fdr != "Pounds" & method.fdr != "Strimmer"
       & method.fdr != "Storey"){
      stop("method.fdr should be one of the 4 methods implemented : BH, Strimmer, Pounds,
           Storey")
    }
    require(fdrtool)
    require(st)
    stud <- studentt.stat(t(exprs(data)),L=t(pData(data)),var.equal=var.equal)
    if(length(stud[is.na(stud)]) != 0){
      nbrNA <- length(stud[is.na(stud)])
      stud[is.na(stud)] <- 0
      warning(nbrNA ," NA's has been replaced by zero")
    }
    
    fdr.stud<-fdrtool(stud,plot=F,verbose=F)
    if(method.fdr=="Pounds"){
      rfdr <- robust.fdr(p=fdr.stud$pval,sides=1,discrete=F)
      names(rfdr$p)<-featureNames(data)
      names.sig.spot<- names(rfdr$p[rfdr$fdr <fdr.thr])
      if(plot){
        plot(rfdr$p,rfdr$fdr, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(method.fdr=="Strimmer"){
      names(fdr.stud$pval)<-featureNames(data)
      names.sig.spot<- names(fdr.stud$pval[fdr.stud$qval <fdr.thr])
      if(plot){
        plot(fdr.stud$pval,fdr.stud$qval, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(method.fdr=="BH"){
      bh.pval<- p.adjust(fdr.stud$pval,method="BH")
      names(bh.pval)<-featureNames(data)
      names.sig.spot<- names(bh.pval[bh.pval <fdr.thr])
      if(plot){
        plot(fdr.stud$pval,bh.pval, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(method.fdr=="Storey"){
      stor<- qvalue(fdr.stud$pval,pi0.method="bootstrap")
      names(stor$qvalues) <- featureNames(data)
      names.sig.spot<- names(stor$qvalues[stor$qvalues <fdr.thr])
      if(plot){
        plot(fdr.stud$pval, stor$qvalues, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(length(names.sig.spot)==0){cat("No significant spots founds")
                                  return(data[0])}
    else{
      a<- data[names.sig.spot]
      a2 <- a[abs(fData(a))>1]
      if(Fold2){
        if(length(exprs(a2))==0){cat("No significant spots founds")
                                 return(data[0])}
        else {
          cat("Number of up-regulated spots in Condition 2\n")
          print(length(featureNames(a2[fData(a2)>0])))
          cat("Number of down-regulated spots in Condition 2\n")
          print(length(featureNames(a2[fData(a2)<0])))
          
          return(a2)}
      }
      else{
        if(length(exprs(a))==0){cat("No significant spots founds")
                                return(data[0])}
        else {
          cat("Number of up-regulated spots in Condition 2\n")
          print(length(featureNames(a[fData(a)>0])))
          cat("Number of down-regulated spots in Condition 2\n")
          print(length(featureNames(a[fData(a)<0])))
          
          return(a)}
      } 
    }
    
    }

#### Smyth's Moderate T-test#####
modT.Prot <-
  function(data, plot = T, fdr.thr = 0.1, Fold2 = F, method.fdr = "BH", col=1){
    if(class(data)[1]!="ExpressionSet"){
      warning(gettextf("'%s' is not an ExpressionSet",
                       as.character(match.call()$data),domain=NA))  
    }
    if(method.fdr != "BH" & method.fdr != "Pounds" & method.fdr != "Strimmer"
       & method.fdr != "Storey"){
      stop("method.fdr should be one of the 4 methods implemented : BH, Strimmer, Pounds,
           Storey")
    }
    
    require(fdrtool)
    require(limma)
    L=t(pData(data)[,col])
    L=as.integer(as.factor(L))
    des =cbind(rep(1, length(L)), L)
    fit <- lmFit(data, design=des,method="robust",maxit=1000)
    fit <- eBayes(fit, robust=T)
    
    if(method.fdr=="Pounds"){
      rfdr <- robust.fdr(p=fit$p.value[,2],sides=1,discrete=T)
      names(rfdr$p)<-featureNames(data)
      names.sig.spot<- names(rfdr$p[rfdr$fdr <fdr.thr])
      if(plot){
        plot(rfdr$p,rfdr$fdr, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(method.fdr=="Strimmer"){
      fdr.lim <- fdrtool(fit$p.value[,2], verbose=F, plot=F, statistic="pvalue")
      names(fdr.lim$pval)<-featureNames(data)
      names.sig.spot<- names(fdr.lim$pval[fdr.lim$qval <fdr.thr])
      if(plot){
        plot(fit$p.value[,2],fdr.lim$lfdr, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(method.fdr=="BH"){
      bh.pval<- p.adjust(fit$p.value[,2],method="BH")
      names(bh.pval)<-featureNames(data)
      names.sig.spot<- names(bh.pval[bh.pval <fdr.thr])
      if(plot){
        plot(fit$p.value[,2],bh.pval, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(method.fdr=="Storey"){
      stor<- qvalue(fit$p.value[,2],pi0.method="bootstrap")
      names(stor$qvalues) <- featureNames(data)
      names.sig.spot<- names(stor$qvalues[stor$qvalues <fdr.thr])
      if(plot){
        plot(fit$p.value[,2], stor$qvalues, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(length(names.sig.spot)==0){cat("No significant spots founds")
                                  return(data[0])}
    else{
      a<- data[names.sig.spot]
      r.2 <- apply(exprs(a[,pData(data)[,col]==levels(pData(data)[,col])[2]]),1,mean)
      r.1 <- apply(exprs(a[,pData(data)[,col]==levels(pData(data)[,col])[1]]),1,mean)
      fData(a)[,1] <- r.2 - r.1
      a2 <- a[abs(fData(a))>1]
      if(Fold2){
        if(length(exprs(a2))==0){cat("No significant spots founds")
                                 return(data[0])}
        else {
          cat("Number of up-regulated spots in Condition 2\n")
          print(length(featureNames(a2[fData(a2)>0])))
          cat("Number of down-regulated spots in Condition 2\n")
          print(length(featureNames(a2[fData(a2)<0])))
          
          return(a2)}
      }
      else{
        if(length(exprs(a))==0){cat("No significant spots founds")
                                return(data[0])}
        else {
          cat("Number of up-regulated spots in Condition 2\n")
          print(length(featureNames(a[fData(a)>0])))
          cat("Number of down-regulated spots in Condition 2\n")
          print(length(featureNames(a[fData(a)<0])))
          
          return(a)}
      } 
    }
    
    }

#### Tusher's modified T-test (same as in Significance Analysis of Microarray)####
samT.Prot <-
  function(data, plot=T,fdr.thr=0.1, Fold2=F, method.fdr="BH"){
    if(class(data)[1]!="ExpressionSet"){
      warning(gettextf("'%s' is not an ExpressionSet",
                       as.character(match.call()$data),domain=NA))  
    }
    if(method.fdr != "BH" & method.fdr != "Pounds" & method.fdr != "Strimmer"
       & method.fdr != "Storey"){
      stop("method.fdr should be one of the 4 methods implemented : BH, Strimmer, Pounds,
           Storey")
    }
    
    require(fdrtool)
    require(samr)
    X=t(exprs(data))
    L=t(pData(data))
    L=as.factor(L)
    dd = list(x = t(X), y = as.integer(L), logged2 = TRUE)
    out =samr(dd, resp.type = "Two class unpaired", nperms = 1,s0=0.9,center.arrays=F)
    fdr.sam <- fdrtool(out$tt, verbose=F, plot=F)
    if(method.fdr=="Pounds"){
      rfdr <- robust.fdr(p=fdr.sam$pval,sides=1,discrete=T)
      names(rfdr$p)<-featureNames(data)
      names.sig.spot<- names(rfdr$p[rfdr$fdr <fdr.thr])
      if(plot){
        plot(rfdr$p,rfdr$fdr, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(method.fdr=="Strimmer"){
      names(fdr.sam$pval)<-featureNames(data)
      names.sig.spot<- names(fdr.sam$pval[fdr.sam$qval <fdr.thr])
      if(plot){
        plot(fdr.sam$pval,fdr.sam$qval, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(method.fdr=="BH"){
      bh.pval<- p.adjust(fdr.sam$pval,method="BH")
      names(bh.pval)<-featureNames(data)
      names.sig.spot<- names(bh.pval[bh.pval <fdr.thr])
      if(plot){
        plot(fdr.sam$pval,bh.pval, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(method.fdr=="Storey"){
      stor<- qvalue(fdr.sam$pval,pi0.method="bootstrap")
      names(stor$qvalues) <- featureNames(data)
      names.sig.spot<- names(stor$qvalues[stor$qvalues <fdr.thr])
      if(plot){
        plot(fdr.sam$pval, stor$qvalues, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(length(names.sig.spot)==0){cat("No significant spots founds")
                                  return(data[0])}
    else{
      a<- data[names.sig.spot]
      a2 <- a[abs(fData(a))>1]
      if(Fold2){
        if(length(exprs(a2))==0){cat("No significant spots founds")
                                 return(data[0])}
        else {
          cat("Number of up-regulated spots in Condition 2\n")
          print(length(featureNames(a2[fData(a2)>0])))
          cat("Number of down-regulated spots in Condition 2\n")
          print(length(featureNames(a2[fData(a2)<0])))
          
          return(a2)}
      }
      else{
        if(length(exprs(a))==0){cat("No significant spots founds")
                                return(data[0])}
        else {
          cat("Number of up-regulated spots in Condition 2\n")
          print(length(featureNames(a[fData(a)>0])))
          cat("Number of down-regulated spots in Condition 2\n")
          print(length(featureNames(a[fData(a)<0])))
          
          return(a)}
      } 
    }
    
    }

### Efron's modified T-test ######
efronT.Prot <-
  function(data, plot=T,fdr.thr=0.1, Fold2=F,method.fdr="BH"){
    if(class(data)[1]!="ExpressionSet"){
      warning(gettextf("'%s' is not an ExpressionSet",
                       as.character(match.call()$data),domain=NA))  
    }
    if(method.fdr != "BH" & method.fdr != "Pounds" & method.fdr != "Strimmer"
       & method.fdr != "Storey"){
      stop("method.fdr should be one of the 4 methods implemented : BH, Strimmer, Pounds,
           Storey")
    }
    require(fdrtool)
    require(st)
    ef<- efront.stat(t(exprs(data)),L=t(pData(data)), verbose=F)
    fdr.ef<-fdrtool(ef,plot=F,verbose=F)
    if(method.fdr=="Pounds"){
      rfdr <- robust.fdr(p=fdr.ef$pval,sides=1,discrete=T)
      names(rfdr$p)<-featureNames(data)
      names.sig.spot<- names(rfdr$p[rfdr$fdr <fdr.thr])
      if(plot){
        plot(rfdr$p,rfdr$fdr, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(method.fdr=="Strimmer"){
      names(fdr.ef$pval)<-featureNames(data)
      names.sig.spot<- names(fdr.ef$pval[fdr.ef$qval <fdr.thr])
      if(plot){
        plot(fdr.ef$pval,fdr.ef$qval, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(method.fdr=="BH"){
      bh.pval<- p.adjust(fdr.ef$pval,method="BH")
      names(bh.pval)<-featureNames(data)
      names.sig.spot<- names(bh.pval[bh.pval <fdr.thr])
      if(plot){
        plot(fdr.ef$pval,bh.pval, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(method.fdr=="Storey"){
      stor<- qvalue(fdr.ef$pval,pi0.method="bootstrap")
      names(stor$qvalues) <- featureNames(data)
      names.sig.spot<- names(stor$qvalues[stor$qvalues <fdr.thr])
      if(plot){
        plot(fdr.ef$pval, stor$qvalues, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(length(names.sig.spot)==0){cat("No significant spots founds")
                                  return(data[0])}
    else{
      a<- data[names.sig.spot]
      a2 <- a[abs(fData(a))>1]
      if(Fold2){
        if(length(exprs(a2))==0){cat("No significant spots founds")
                                 return(data[0])}
        else {
          cat("Number of up-regulated spots in Condition 2\n")
          print(length(featureNames(a2[fData(a2)>0])))
          cat("Number of down-regulated spots in Condition 2\n")
          print(length(featureNames(a2[fData(a2)<0])))
          
          return(a2)}
      }
      else{
        if(length(exprs(a))==0){cat("No significant spots founds")
                                return(data[0])}
        else {
          cat("Number of up-regulated spots in Condition 2\n")
          print(length(featureNames(a[fData(a)>0])))
          cat("Number of down-regulated spots in Condition 2\n")
          print(length(featureNames(a[fData(a)<0])))
          
          return(a)}
      } 
    }
    
    }


### Strimmer's Shrinked T-test ######
shrinkT.Prot <-
  function(data, plot = T, fdr.thr = 0.1, Fold2 = F, method.fdr = "BH", var.equal = F){
    if(class(data)[1]!="ExpressionSet"){
      warning(gettextf("'%s' is not an ExpressionSet",
                       as.character(match.call()$data),domain=NA))  
    }
    if(method.fdr != "BH" & method.fdr != "Pounds" & method.fdr != "Strimmer"
       & method.fdr != "Storey"){
      stop("method.fdr should be one of the 4 methods implemented : BH, Strimmer, Pounds,
           Storey")
    }
    require(fdrtool)
    require(st)
    shrt <- shrinkt.stat(t(exprs(data)),L=t(pData(data)), 
                         var.equal=var.equal,verbose=F)
    fdr.shrt<-fdrtool(shrt,plot=F,verbose=F)
    if(method.fdr=="Pounds"){
      rfdr <- robust.fdr(p=fdr.shrt$pval,sides=1,discrete=T)
      names(rfdr$p)<-featureNames(data)
      names.sig.spot<- names(rfdr$p[rfdr$fdr <fdr.thr])
      if(plot){
        plot(rfdr$p,rfdr$fdr, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(method.fdr=="Strimmer"){
      names(fdr.shrt$pval)<-featureNames(data)
      names.sig.spot<- names(fdr.shrt$pval[fdr.shrt$qval <fdr.thr])
      if(plot){
        plot(fdr.shrt$pval,fdr.shrt$qval, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(method.fdr=="BH"){
      bh.pval<- p.adjust(fdr.shrt$pval,method="BH")
      names(bh.pval)<-featureNames(data)
      names.sig.spot<- names(bh.pval[bh.pval <fdr.thr])
      if(plot){
        plot(fdr.shrt$pval,bh.pval, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(method.fdr=="Storey"){
      stor<- qvalue(fdr.shrt$pval,pi0.method="bootstrap")
      names(stor$qvalues) <- featureNames(data)
      names.sig.spot<- names(stor$qvalues[stor$qvalues <fdr.thr])
      if(plot){
        plot(fdr.shrt$pval, stor$qvalues, xlab="p-value",ylab="Fdr")
        abline(h=fdr.thr)
      }
    }
    if(length(names.sig.spot)==0){cat("No significant spots founds")
                                  return(data[0])}
    else{
      a<- data[names.sig.spot]
      a2 <- a[abs(fData(a))>1]
      if(Fold2){
        if(length(exprs(a2))==0){cat("No significant spots founds")
                                 return(data[0])}
        else {
          cat("Number of up-regulated spots in Condition 2\n")
          print(length(featureNames(a2[fData(a2)>0])))
          cat("Number of down-regulated spots in Condition 2\n")
          print(length(featureNames(a2[fData(a2)<0])))
          
          return(a2)}
      }
      else{
        if(length(exprs(a))==0){cat("No significant spots founds")
                                return(data[0])}
        else {
          cat("Number of up-regulated spots in Condition 2\n")
          print(length(featureNames(a[fData(a)>0])))
          cat("Number of down-regulated spots in Condition 2\n")
          print(length(featureNames(a[fData(a)<0])))
          
          return(a)}
      } 
    }
    
    }


###### Simulation of 2D Gel Volume data #####
Sim.Prot.2D <- 
  function (data,nsp=nrow(data),nr=10,p0=0.1, s2_0=0.2,d0=3){
    if(class(data)!="data.frame" & class(data)!="matrix"){
      warning(gettextf("'%s' is neither a dataframe nor a matrix",
                       as.character(match.call()$data),domain=NA))
    }
    if(nsp > nrow(data)){
      stop(gettextf("nsp must be equal or less than the number of rows in '%s'",
                    as.character(match.call()$data),domain=NA))
    }
    
    Cond <- data.frame(Cond=c(rep("Cond1",nr),rep("Cond2",nr)))
    rownames(Cond) <- c(as.vector(paste("A",as.character(1:nr),sep="")),
                        as.vector(paste("B",as.character(1:nr),sep="")))  
    
    s2_g=d0 * s2_0/rchisq(nsp, df = d0)
    dat=matrix(0,nrow=nsp,ncol=nr+nr)
    rownames(dat)<- as.character(c(1:nsp))
    ndx1 = runif(nsp) < p0/2
    nde1 = length(which(ndx1))
    nam1= which(ndx1)
    ndx2 = runif(nsp) < p0/2
    nde2 = length(which(ndx2)) 
    nam2= which(ndx2)
    
    mean.in <- apply(data, 1, mean)
    
    for(n in 1:nsp){
      dat[n,] = matrix(rnorm(nr+nr, mean = mean.in[n], sd =sqrt(s2_g)))  }
    
    for(n in 1:nde1){
      dat[nam1[n],1:nr]= matrix(rnorm(nr, mean = mean.in[nam1[n]]+1, 
                                      sd = 2*sqrt(s2_g)))  }
    for(n in 1:nde2){
      dat[nam2[n],c(nr+1):c(nr+nr)]= matrix(rnorm(nr, mean = mean.in[nam2[n]]+1,
                                                  sd = 2*sqrt(s2_g)))  }
    Sim<-dat
    colnames(Sim) <- rownames(Cond)
    ES.sim <- ES.prot(Sim,n1=nr,n2=nr,f=Cond)
    notes(ES.sim)$SpotSig <- as.character(c(nam1,nam2))
    
    return(ES.sim)
  }

