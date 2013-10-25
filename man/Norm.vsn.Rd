\name{Norm.vsn}
\alias{Norm.vsn}
\title{
Normalize 2D Gel Volume data using VSN
}
\description{
This function allows to normalize 2D Gel Volume data using the "Variance Stabilizing Normalization".
}
\usage{
Norm.vsn(data, n1, n2, plot = T)
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
The "Variance Stabilizing Normalization" relies on a transformation \emph{h}, of the parametric form \emph{h(x)= arsinh(a+bx)}  (for details see, Huber et al., 2002). The parameters of \emph{h} together with those of the calibration between experiments are estimated with a robust variant of maximum-likelihood estimation.This function is based on \code{\link{normalizeVSN}} from \code{limma} package.
}

\value{
The function returns a matrix of vsn normalized data
}
\author{
Sebastien Artigaud
\email{sebastien.artigaud@gmx.com}
}
\references{
\itemize{
\item{Artigaud, S., Gauthier, O. & Pichereau, V. (2013) "Identifying differentially expressed proteins in two-dimensional electrophoresis experiments: inputs from transcriptomics statistical tools." Bioinformatics, vol.29 (21): 2729-2734.}
\item{
Huber, W., Heydebreck, von, A., Sultmann, H., Poustka, A., & Vingron, M. (2002) "Variance stabilization applied to microarray data calibration and to the quantification of differential expression" Bioinformatics, vol. 18 (Suppl 1): S96-S104.}
}
}
\seealso{
\code{\link{normalizeVSN}},\code{\link{RIplot}},\code{\link{ES.prot}}
}
\examples{
data(pecten)

pecten.norm <- Norm.vsn(pecten, n1=6, n2=6, plot=TRUE)
}
\keyword{prot2D}
