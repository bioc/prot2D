\name{Sim.Prot.2D}
\alias{Sim.Prot.2D}
\title{
Simulate 2D Gel volume data
}
\description{
This function simulates realistic 2D Gel volume data, based on parameters estimate from real dataset.
}
\usage{
Sim.Prot.2D(data, nsp=nrow(data), nr=10, p0=0.1, s2_0=0.2 ,d0=3)
}
\arguments{
  \item{data}{
either a dataframe or a matrix of 2D Gel Volume data. Usually the output of \code{\link{Norm.qt}}.
}
  \item{nsp}{
an integer. Number of spots (i.e. rows) to generate (should be equal or less than the number of rows of data).
}
  \item{nr}{
an integer. Number of replicates for each conditions.
}
  \item{p0}{
proportion of spots to be differentially expressed (defaults to 0.1).
}
  \item{s2_0}{
estimator of the standard deviation of the spots distribution (defaults to 0.2).
}
  \item{d0}{
degree of freedom of the Chi-squared distribution used in the calculation of the distribution (defaults to 3).
}
}
\details{
Volume data are computed following these steps (see Smyth, 2004 and Artigaud et al, 2013 for details) :
\enumerate{
\item{Log2 mean volumes from data are computed for each spot.}
\item{Means are used as input parameters in order to simulate a normal distribution (with no differential expression between conditions) for each spots with standard deviations, computed as described by Smyth (2004).}
\item{p0 of the spots are randomly picked for introducing differential expression in both conditions (p0/2 in each condition).}
}
If \code{nsp} is less than the number of rows in \code{data}, the \code{nsp} first columns of data are used as input parameters.
}
\value{
returns an \code{\link{ExpressionSet}} of simulated volume data (log2 transformed) with 2 conditions ("Cond1" and "Cond2" in \code{\link{phenoData}}) slot of the \code{\link{ExpressionSet}}.The spots differentially generated can be retrieve with \code{\link{notes}}.
}
\references{
\itemize{
\item{Smyth, G.K. (2004) "Linear models and empirical bayes methods for assessing differential expression in microarray experiments." Statistical Applications in Genetics and Molecular Biology, vol. 3: Article 3.}
\item{Artigaud, S., Gauthier, O. & Pichereau, V. (2013) "Identifying differentially expressed proteins in 2-DE experiments: Inputs from transcriptomics statistical tools." Submitted}
}
}

\author{
Sebastien Artigaud
\email{sebastien.artigaud@gmx.com}
}
\seealso{
\code{\link{variousT.Prot}}, 
}
\examples{
#Simulate data based on "pecten"
data(pecten)
data(pecten.fac)

pecten.norm <- Norm.qt(pecten, n1=6, n2=6, plot=FALSE)
Sim.data <- Sim.Prot.2D(data=pecten.norm, nsp=700, nr=10, p0=0.1, s2_0=0.2, d0=3)

#Compare different methods for finding diffenrentially expressed proteins
res.stud <- ttest.Prot(Sim.data, fdr.thr=0.1, plot=FALSE)
res.mo <- modT.Prot(Sim.data, fdr.thr=0.1, plot=FALSE)

a <- featureNames(res.stud) 
a # Names of the spots selected for a student's t-test with an FDR of 0.1
b <- featureNames(res.mo) 
b # Names of the spots selected for a  modT-test with an FDR of 0.1
c <- notes(Sim.data)$SpotSig 
c #  Names of the really differentially generated spots

# Compare between t-test and moderate t-test (requires made4 package)
# Just for illustration purpose (only 1 simulation !)
\dontrun{
cl.t <- comparelists(c,a)
cl.m <- comparelists(c,b)

fp <- matrix(c(c(length(a)- (cl.t$Length.Y-length(cl.t$intersect))),
c(cl.t$Length.Y-length(cl.t$intersect)), 
c(length(b)- (cl.m$Length.Y-length(cl.m$intersect))),
c(cl.m$Length.Y-length(cl.m$intersect))),
nrow=2,ncol=2,byrow=F)

barplot(fp, main="spots declared significant",names.arg=c("t-test","modT"))
legend("topleft",c("False Positives","True Positives"), fill=c("grey","black"), bty="n")
}

} 

\keyword{prot2D}
