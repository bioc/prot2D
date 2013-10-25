\name{prot2D-package}
\alias{prot2D-package}
\alias{prot2D}
\docType{package}
\title{Statistics Tools for data issued from a 2D Gel Electrophoresis}
\description{
 The purpose of this package is to analyze (i.e. Normalize 
and select significant spots) data issued from 2D GEl experiments
}
\details{
This package provides a simple interface for analysing data from 2D gel experiments. Functions for performing normalization as well as selecting significant spots are provided. All the functions for selecting significant spots are adapted from functions for microarray analysis provided by Bioconductor and CRAN, and all credits go to the authors of these functions. For analyzing 2D gel experiments data, users are advised to follow theses steps :
\enumerate{
\item{Normalize data using \code{\link[prot2D]{Norm.qt}}}
\item{Coerce data into an  \code{\link[Biobase]{ExpressionSet}} using \code{\link[prot2D]{ES.prot}}}
\item{Use normalized data to find differentially expressed proteins with FDR-controled functions: \code{\link[prot2D]{modT.Prot}} was find to be the more efficient for 2-DE (see Artigaud et al , 2013) but other functions are provided (\code{\link{ttest.Prot}},\code{\link{samT.Prot}},\code{\link{efronT.Prot}}, ,\code{\link{shrinkT.Prot}}) }
}
}
\author{
Sebastien Artigaud
\email{sebastien.artigaud@gmx.com}
}
\references{
\itemize{
\item{Artigaud, S., Gauthier, O. & Pichereau, V. (2013) "Identifying differentially expressed proteins in two-dimensional electrophoresis experiments: inputs from transcriptomics statistical tools." Bioinformatics, vol.29 (21): 2729-2734.}
\item{Dudoit, S., Yang, Y.H., Callow, M.J., & Speed, T.P. (2002) "Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments" Statistica Sinica, vol. 12: 111-139.}
\item{Strimmer, K. (2008) "A unified approach to false discovery rate estimation." BMC Bioinformatics, vol. 9: 303.}
}
\seealso{
\code{\link{ExpressionSet}},
\code{\link{fdrtool}}
}
\examples{
data(pecten)
data(pecten.fac)

pecten.norm <- Norm.qt(pecten, n1=6, n2=6, plot=TRUE) #Quantiles normalization of the data
ES.p <- ES.prot(pecten.norm, n1=6, n2=6, f=pecten.fac)
x <- modT.Prot(ES.p, fdr.thr=0.1, plot=TRUE)
featureNames(x) # Names of the spots selected for a moderated t-test with a fdr of 0.1
fData(x) # Displaying fold change (as log2(ratio)) for selected spots
exprs(x) # Normalized volume data for all the selected spots
\dontrun{heatplot(x)} #Great heatmap of the selected spots (requires made4 Bioconductor package )
}
\keyword{ package }
\keyword{ biology }
\keyword{ proteomic }
\keyword{ 2DGel }
\keyword{prot2D}
