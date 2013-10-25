\name{Norm.qt}
\alias{Norm.qt}
\title{
Normalize 2D Gel Volume data using Quantiles Normalization
}
\description{
This function allows to normalize 2D Gel Volume data using the Quantiles Normalization.
}
\usage{
Norm.qt(data, n1, n2, plot = T)
}
\arguments{
  \item{data}{
a dataframe of raw 2D Gel Volume data. \code{data} should be raw intensities displayed with gel as columns with the name of columns corresponding to the names of the gels and spots as rows with the names of the rows corresponding to the name of the spots. The replicates for each condition should be ordered in following columns.
}
  \item{n1}{
an integer. Number of replicates in condition 1.
}
  \item{n2}{
an integer. Number of replicates in condition 2.
}
  \item{plot}{
logical. if TRUE (default) displaying two \code{\link{RIplot}}, one with the raw data, another with normalized data.
}
}
\details{
2D Gel Volume data must be normalized in order to remove systemic variation prior to data analysis. The principle of the "quantiles normalization" is to set each quantile of each column (i.e. the spots volume data of each gels)  to the mean of that quantile across gels. The intention is to make all the normalized columns have the same empirical distribution.This function is based on \code{\link{normalizeQuantiles}} from \code{limma} package.
}

\value{
The function returns a matrix of log2 transformed quantiles normalized data
}
\author{
Sebastien Artigaud
\email{sebastien.artigaud@gmx.com}
}
\references{
Artigaud, S., Gauthier, O. & Pichereau, V. (2013) "Identifying differentially expressed proteins in two-dimensional electrophoresis experiments: inputs from transcriptomics statistical tools." Bioinformatics, vol.29 (21): 2729-2734.
}
\seealso{
\code{\link{normalizeQuantiles}},\code{\link{RIplot}},\code{\link{ES.prot}}
}
\examples{
data(pecten)

pecten.norm <- Norm.qt(pecten, n1=6, n2=6, plot=TRUE)
}
\keyword{prot2D}
