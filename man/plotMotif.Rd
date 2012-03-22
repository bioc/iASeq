\name{plotMotif}
\alias{plotMotif}
\title{Correlation Motif Plot}
\description{This function plots the Correlation Motif patterns, the associated prior probability distributions and the number of SNPs called for each motif based on posterior probability.}

\usage{
plotMotif(bestmotif,title="",cutoff)
}

\arguments{
\item{bestmotif}{	The bestmotif obtained from iASeqmotif. }
\item{title}{		The title on the figure. }
\item{cutoff}{	The posterior probability cutoff for calling a SNP belonging to certain motif. }
}

\details{Each row in all graphs corresponds to one motif pattern. The first graph shows \eqn{qup}, the correlation motif pattern of over expression (binding). The second graph shows \eqn{qdown}, the correlation motif pattern of under expression (binding). The grey color scale of cell \eqn{(k,d)} indicates the probability that motif \eqn{k} is over or under expressed in study \eqn{d}. Each row of the two bar charts corresponds
to the motif pattern in the same row of the left two pattern graphs. The length of
the bar in the first bar chart estimates the number of SNPs of the given pattern in the
dataset according to motif frequency, which is equal to \eqn{motif.fitted\$bestmotif\$motif.prior} multiplying
the number of total SNPs. The length of the bar in the second bar chart shows the number of SNPs called for the given pattern according to the \eqn{cutoff} of posterior probability.}

\author{Yingying Wei, Hongkai Ji}

\references{Yingying Wei, Xia Li, Qianfei Wang, Hongkai Ji (2012) iASeq: integrating multiple ChIP-seq datasets for detecting allele-specific binding.}


\seealso{
\code{\link{iASeqmotif}}, \code{\link{plotBIC}}, \code{\link{sampleASE}}
}

\examples{
example(iASeqmotif) # compute 'motif.fitted'
plotMotif(motif.fitted$bestmotif,cutoff=0.9)
}


