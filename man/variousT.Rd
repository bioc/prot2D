\name{variousT.Prot}
\alias{variousT.Prot}
\alias{ttest.Prot}
\alias{modT.Prot}
\alias{samT.Prot}
\alias{efronT.Prot}
\alias{shrinkT.Prot}
\title{
Various Tests for finding differentially expressed proteins in 2-DE experiments
}
\description{
These functions provides simple methods for finding differentially expressed proteins in 2-DE experiments.
}

\usage{
ttest.Prot(data, plot = T, fdr.thr = 0.1, Fold2 = F, method.fdr = "BH", var.equal = F)

modT.Prot(data, plot = T, fdr.thr = 0.1, Fold2 = F, method.fdr = "BH", col=1)

samT.Prot(data, plot = T, fdr.thr = 0.1, Fold2 = F, method.fdr = "BH")

efronT.Prot(data, plot = T, fdr.thr = 0.1, Fold2 = F, method.fdr = "BH")

shrinkT.Prot(data, plot = T, fdr.thr = 0.1, Fold2 = F, method.fdr = "BH", var.equal = F)
}
\arguments{
  \item{data}{ an \code{\link{ExpressionSet}} of volume data. Usually as returned by \code{\link{ES.prot}}.
}
  \item{fdr.thr}{threshold value for the FDR (defaults to 0.1). 
}
  \item{method.fdr}{method for FDR estimator. 4 methods are implemented :"BH" for the Classical FDR from Benjamini & Hochberg, "Strimmer", "Storey" and "Pounds" (see references for details).
}
  \item{Fold2}{logical. If true only spots with an absolute ratio of 2 are selected.
}
  \item{plot}{logical. If true (default) draws a plot fdr values depending on the p-values with a line indicating the FDR threshold. 
}
 \item{var.equal}{logical. If false (default) assume unequal variances in each condition and uses a proper correction. Only for \code{ttest.Prot} and \code{shrinkT.Prot}.
 }
   \item{col}{An integer (default to 1). Only for \code{modT.Prot}. If there is more than one factor in \code{\link{pData}}, indicates the column to use for the analysis}.
}
\details{
As described in Artigaud et al (2013) these are functions adapted from microarray analysis. 2-DE experiments analysis requires a variant of the t-statistic that is suitable for high-dimensional data and large-scale multiple testing. For this purpose, in the last few years, various test procedures have been suggested.\cr
These functions provides:\cr
 - the classical Student's t-test (adapted from \code{\link{studentt.stat}}).\cr
 - two tests especially modified for micro-array analysis : Efron's t-test (adapted from \code{\link{efront.stat}}, Efron et al, 2001) and the modified t-test used in Significance Analysis for Microarray (adapted from \code{\link{samr}}, Tusher et al, 2001)\cr
 - two methods that take advantage of hierarchical Bayes methods for estimation of the variance across genes: the moderate t-test from Smyth (using \code{\link{limma}}; see Smyth, 2004) and the "Shrinkage t" statistic test from Opgen-Rhein & Strimmer (adapted from \code{\link{shrinkcat.stat}}; see Opgen-Rhein & Strimmer, 2007). \cr
As statistical tests allowing the identification of differentially expressed proteins must take into account a correction for multiple tests in order to avoid false conclusions. These functions also provides different methods to estimate the False Discovery Rate :\cr
- the classical FDR estimator of Benjamini & Hochberg (using \code{\link{p.adjust}}; see Benjamini & Hochberg, 1995)\cr
- the Fdr estimator of Strimmer (based on local fdr calculation) (using \code{\link{fdrtool}}; see Strimmer 2008)\cr
- the "robust FDR" estimator of Pounds & Cheng (implemented in \code{\link{robust.fdr}} for the \code{\link{prot2D}} package; see Pounds & Cheng, 2006)\cr
- Fdr method of Storey and Tibshirani (2003), also known as "q-values" (using \code{\link{qvalue}}.
}
\value{
returns an \code{\link{ExpressionSet}} containing only the significant spots (see Examples).

}
\references{
\itemize{
\item{Artigaud, S., Gauthier, O. & Pichereau, V. (2013) "Identifying differentially expressed proteins in two-dimensional electrophoresis experiments: inputs from transcriptomics statistical tools." Bioinformatics, vol.29 (21): 2729-2734.}
\item{Benjamini, Y. & Hochberg, Y. (1995) "Controlling the false discovery rate: a practical and powerful approach to multiple testing" Journal of the Royal Statistical Society. Series B. Methodological.: 289-300.}
\item{Efron, B., Tibshirani, R., Storey, J.D., & Tusher, V. (2001) "Empirical Bayes Analysis of a Microarray Experiment" Journal of the American Statistical Association, vol. 96 (456): 1151-1160.}
\item{Tusher, V.G., Tibshirani, R., & Chu, G. (2001) "Significance analysis of microarrays applied to the ionizing radiation response"" Proceedings of the National Academy of Sciences of the United States of America, vol. 98 (9): 5116-5121.}
\item{Smyth, G.K. (2004) "Linear models and empirical bayes methods for assessing differential expression in microarray experiments." Statistical Applications in Genetics and Molecular Biology, vol. 3: Article 3.}
\item{Pounds, S. & Cheng, C. (2006) "Robust estimation of the false discovery rate" Bioinformatics, vol. 22 (16): 1979-1987.}
\item{Strimmer, K. (2008) "A unified approach to false discovery rate estimation." BMC Bioinformatics, vol. 9: 303.}
\item{Opgen-Rhein, R. & Strimmer, K. (2007) "Accurate Ranking of Differentially Expressed
Genes by a Distribution-Free Shrinkage Approach" Statistical Applications in Genetics and Molecular Biology, vol. 6 (1).}
}

\author{
Sebastien Artigaud
\email{sebastien.artigaud@gmx.com}
}

\seealso{
\code{\link{Norm.qt}},\code{\link{ES.prot}},\code{\link{fdrtool}},\code{\link{limma}},
\code{\link{samr}},\code{\link{studentt.stat}},\code{\link{shrinkt.stat}},
\code{\link{efront.stat}},\code{\link{qvalue}}
}
\examples{
data(pecten)
data(pecten.fac)

pecten.norm <- Norm.qt(pecten, n1=6, n2=6, plot=TRUE) #Quantiles normalization of the data
ES.p <- ES.prot(pecten.norm, n1=6, n2=6, f=pecten.fac)
ES.diff <- modT.Prot(ES.p, fdr.thr=0.1, plot=TRUE)
featureNames(ES.diff) # Names of the spots selected for a moderated t-test with a fdr of 0.1
fData(ES.diff) # Displaying fold change (as log2(ratio)) for selected spots
exprs(ES.diff) # Volume normalized data for all the selected spots
\dontrun{heatplot(ES.diff)} #Great heatmap of the selected spots (require made4 Bioconductor package )}
\keyword{prot2D}
