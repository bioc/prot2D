\name{Norm.loess}
\alias{Norm.loess}
\title{
Normalize 2D Gel Volume data using loess
}
\description{
This function allows to normalize 2D Gel Volume data between two conditions using
\code{\link{loess}} smoothing.
}
\usage{
Norm.loess(data, n1, n2, plot = T, span = 0.75)
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
logical. if TRUE (default) displaying two \code{\link{RIplot}}, one with the raw data, another with \code{\link{loess}} normalized data.
}
  \item{span}{
parameter for \code{\link{loess}} regression, which sets the proportion of the total data set to be used in each window for local fitting (defaults to 0.75)
}
}
\details{
2D Gel Volume data must be normalized in order to remove systemic variation prior to data analysis (for details, see Artigaud et al, 2013). The loess method is a robust smoother based on local polynomial regression, raw spot's volumes data are corrected by
\deqn{V'_{Cond1} = \frac{V_{Cond1}}{\sqrt{2^{-loess fit}}}}{%
V'_{Cond1} = V_{Cond1} / \sqrt(2^- loess fit)}
\deqn{V'_{Cond2} = V_{Cond2}\times\sqrt{2^{-loess fit}}}{%
V'_{Cond2} = V_{Cond2} x \sqrt(2^- loess fit)}
}

\value{
The function returns a matrix of loess normalized data
}
\author{
Sebastien Artigaud
\email{sebastien.artigaud@gmx.com}
}
\references{
Artigaud, S., Gauthier, O. & Pichereau, V. (2013) "Identifying differentially expressed proteins in 2-DE experiments: Inputs from transcriptomics statistical tools." Submitted
}
\seealso{
\code{\link{loess}},\code{\link{RIplot}},\code{\link{ES.prot}}
}
\examples{
data(pecten)
data(pecten.fac)
pecten.norm <- Norm.loess(pecten, n1=6,n2=6,plot=TRUE,span=0.75)
}
\keyword{prot2D}
