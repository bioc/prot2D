\name{robust.fdr}
\alias{robust.fdr}
\title{
"robust FDR" estimation 
}
\description{
Implements robust method of FDR estimation (Pounds and Cheng 2006, Bioinformatics)
}
\usage{
robust.fdr(p, sides = 1, p2 = 1 - p, discrete = F, use8 = T)
}

\arguments{
  \item{p}{vector of p-values from the analysis.}
  \item{sides}{indicate whether p-values are 1-sided (set sides=1) or 2-sided (set sides=2), default=1.}
  \item{p2}{for one-sided testing, p-values from testing the "other alternative", default=1-p.}
  \item{discrete}{logical. Indicates whether p-values are discrete}
  \item{use8}{indicates whether the constant 8 should be used if p-values are discrete, see Pounds and Cheng (2006) for more details.}
}
\details{
This function uses the code from Stan Pounds available at \url{http://www.stjuderesearch.org/depts/biostats/documents/robust-fdr.R} and is included in \code{\link{prot2D}} package for convenience and comparison purpose.
}
\value{
A list with components:
\item{p}{the vector of p-values provided by the user.}
\item{fdr}{the vector of smoothed FDR estimates.}
\item{q}{the vector of q-values based on the smoothed FDR estimates.}
\item{cdf}{the vector with p-value empirical distribution function at corresponding entry of p.}
\item{loc.fdr}{the local (unsmoothed) FDR estimates.}
\item{fp}{the estimated number of false positives at p-value cutoff in p.}
\item{fn}{the estimated number of false negatives at p-value cutoff in p.}
\item{te}{the total of fp and fn.}
\item{pi}{the null proportion estimate.}
\item{ord}{a vector of indices to order the vectors above by ascending p-value.}
}
\references{
Pounds, S. & Cheng, C. (2006) "Robust estimation of the false discovery rate" Bioinformatics, vol. 22 (16): 1979-1987.
}
\author{
Stan Pounds. Edited by Sebastien Artigaud for \code{\link{prot2D}} package.
}

\keyword{prot2D}