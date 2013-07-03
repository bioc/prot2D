\name{ES.prot}
\alias{ES.prot}
\title{
Create an ExpressionSet containing 2D Gel Volume Data
}
\description{
This function create an an \code{\link{ExpressionSet}} containing 2D Gel Volume Data and associated metadata. It requires a matrix of 2D Gel Volume data with Gels as columns and Spots as rows as well as a dataframe describing the data.}
\usage{
ES.prot(data, n1, n2, f)
}
\arguments{
  \item{data}{
a matrix of spots intensities. \code{data} should be intensities displayed with gels as columns with the name of columns corresponding to the names of the gels and spots as rows with the names of the rows corresponding to the name of the spots. The replicates for each condition should be ordered in following columns.
}
  \item{n1}{
an integer. Number of replicates in condition 1.
}
  \item{n2}{
an integer. Number of replicates in condition 2.
}
  \item{f}{
a dataframe giving the factors for data. \code{f} dataframe should have only 1 column giving the 2 levels of factor with \code{rownames} corresponding to the names of gels (in the same order as \code{data}).
}
}
\details{
Usually \code{data} matrix should be the output of \code{\link{Norm.qt}} and therefore is log2-transformed. In \code{data} the columns (i.e. Gels) must be ordered with replicates from condition 1 before those from condition 2.
}
\value{
The function returns an \code{\link{ExpressionSet}}. The matrix of spots intensities (i.e. Volume) is stored in the \code{\link{assayData}} slot of the \code{\link{ExpressionSet}} and can be retrieved with \code{\link{exprs}}. In the \code{\link{featureData}} slot the log2-ratio is computed and can be retrieved with \code{\link{fData}}.
}
\author{
Sebastien Artigaud
\email{sebastien.artigaud@gmx.com}
}
\seealso{
\code{\link{Norm.qt}},\code{\link{Norm.vsn}},\code{\link{ExpressionSet}}

}
\examples{
data(pecten)
data(pecten.fac)

pecten.norm <- Norm.qt(pecten, n1=6, n2=6, plot=TRUE) #Quantiles normalization of the data
ES.p <- ES.prot(pecten.norm, n1=6, n2=6, f=pecten.fac)
head(exprs(ES.p))}
\keyword{prot2D}
