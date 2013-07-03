\name{RIplot}
\alias{RIplot}
\title{
Ratio-Intensity Plot for 2D Gel Volume data
}
\description{
R-I plots consist in plotting the intensity log2-ratio (R) against mean log10 intensity (I) of spot's volume data from 2D Gel experiments.
}
\usage{
RIplot(data, n1, n2, ...)
}
\arguments{
  \item{data}{
either a dataframe or a matrix of 2D Gel Volume data. \code{data} should be spots intensities displayed with gel as columns with the name of columns corresponding to the names of the gels and spots as rows with the names of the rows corresponding to the name of the spots. The replicates for each condition should be ordered in following columns.
}
  \item{n1}{
an integer. Number of replicates in condition 1.
}
  \item{n2}{
an integer. Number of replicates in condition 2.
}
  \item{...}{Further arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}).}
}
\details{
Dudoit et al (2002) proposed a method for visualization of artifacts in microarray datasets, called the MA-plot, which was transposed for proteomics data as the Ratio-Intensity plot by Meunier et al (2005). R-I plots allow to directly visualize artifacts in the original data set as well as the effects of normalization. The log2-ratio (R) is the log2 of the mean of volume data in condition 2 upon the mean of volume data in condition 1:
\deqn{R = log_2\frac{mean(V_{Cond2})}{mean(V_{Cond1})}}{%
R = log2(mean(V_{Cond2})/mean(V_{Cond1}))}
The intensity (I) is the log10 of the mean of volume data in condition 2 by the mean of volume data in condition 1
\deqn{I = log_10(mean(V_{Cond2})\times mean(V_{Cond1}))}{%
I = log10(mean(V_{Cond2}) x mean(V_{Cond1}))}
}
\value{
Display the RI-plot of data, the values on the plot are the number of spots with a Ratio greater then 1 (upper values) or lesser then -1 (lower value).
A list is invinsibly returns with components :
\item{RI}{a dataframe with ratio and intensity values calculated.}
\item{diff}{a vector with number of spots with a Ratio greater then 1 or lesser then -1 (in that order).}
}
\references{
\itemize{
\item{Dudoit, S., Yang, Y.H., Callow, M.J., & Speed, T.P. (2002) "Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments" Statistica Sinica, vol. 12: 111-139.}
\item{Meunier, B., Bouley, J., Piec, I., Bernard, C., Picard, B., & Hocquette, J.-F. (2005) "Data analysis methods for detection of differential protein expression in two-dimensional gel electrophoresis" Analytical Biochemistry, vol. 340 (2): 226-230.}
}
}
\author{
Sebastien Artigaud
\email{sebastien.artigaud@gmx.com}
}
\seealso{
\code{\link{Norm.qt}},\code{\link{Norm.vsn}}
}
\examples{
data(pecten)

RIplot(pecten, 6, 6, main="Without Normalization")
}

\keyword{prot2D}