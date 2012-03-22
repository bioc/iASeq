\name{iASeqmotif}
\alias{iASeqmotif}
\title{Correlation Motif Fit for Allele Specific Events}
\description{This function fits the Correlation Motif model to multiple RNAseq or ChIPseq studies. It gives the fitted values for the probability distribution of each motif, the fitted values of the given correlation matrix and the posterior probability for each SNP to be allele-specific events (allele-specific expression or allele-specific binding).}

\usage{
iASeqmotif(exprs,studyid,repid,refid,K,iter.max=100,tol=1e-3)
}

\arguments{
\item{exprs}{ A matrix, each row of the matrix corresponds to a heterozygotic SNP and each column of the matrix corresponds to the reads count for either the reference allele or non-reference allele in a replicate of a study. }
\item{studyid}{ The group label for each column of exprs matrix. all columns in the same study have the same studyid. }
\item{repid}{ The sample label for each column of exprs matrix. The two columns within the same sample, one for reference allele and the other for non-reference allele, have the same repid. In other words, repid discriminates the different replicates within the same study.}
\item{refid}{ The reference allele label for each column of exprs matrix. Please code 0 for reference allele columns and 1 for non-reference allele columns to make the interpretation of over expressed(or bound) to be skewing to the reference allele. Otherwise, just interpret the other way round.}
\item{K}{ A vector, each element specifing the number of non-null motifs a model wants to fit. }
\item{tol}{ The relative tolerance level of error. }
\item{iter.max}{ Maximun number of iterations. }
}

\details{For the i^{th} element of \eqn{K}, the function fits total number of \eqn{K[i]+1} motifs, \eqn{K[i]} non-null motifs and the null motif, to the data. Each SNP can belong to one of the \eqn{K[i]+1} possible motifs according to prior probability distribution, \eqn{motif.prior}. For SNPs in motif \eqn{j} \eqn{(j>=1)}, the probability that they are over expressed (or bound) for the reference allele in study \eqn{d} is \eqn{motif.qup(j,d)} and the probability that they are under expressed (or bound) is \eqn{motif.qdown(j,d)}. One should indicate the studyid, repid and refid for each column clearly.
}


\value{
\item{bestmotif$p.post}{ The posterior probability for each SNP to be allele-specific event. A vector whose length correpsonds to the number of SNPs.}
\item{bestmotif$motif.prior}{ Fitted values of the probability distribution of the \eqn{K[i]+1}  motifs for the best fitted model, the first element specifies the null motif and the 2nd to \eqn{K[i]+1}th element correspond to the \eqn{K[i]} non-null motifs.}
\item{bestmotif$motif.qup}{ Fitted values of the over expressed (or bound) correlation motif matrix for the best fitted model. Each row corresponds to a non-null motif and each column corresponds to a study.}
\item{bestmotif$motif.qdown}{ Fitted values of the under expressed (or bound) correlation motif matrix for the best fitted model. Each row corresponds to a non-null motif and each column corresponds to a study.}
\item{bestmotif$clustlike}{ Posterior probability for a SNP to belong to a specific motif based on the best fitted model. Each row corresponds to a SNP and each column corresponds to a motif class.}
\item{bestmotif$c0j}{ \eqn{\alpha} parameter for the null beta prior distribution for each sample.}
\item{bestmotif$d0j}{ \eqn{\beta} parameter for the null beta prior distribution for each sample.}
\item{bestmotif$loglike}{ The log-likelihood for the best fitted model.}
\item{bic}{ The BIC values of all fitted models. A matrix whose first column is the same as input motif number vector ('K') and the second column corresponds to the BIC value of model given by the motif number in the first column in the same row.}
\item{loglike}{ The log-likelihood of all fitted models. A matrix whose first column is the same as input motif number vector ('K') and the second column corresponds to the log likelihood value of the model given by the motif number in the first column in the same row.}
}

\seealso{
\code{\link{plotBIC}}, \code{\link{plotMotif}}, \code{\link{sampleASE}}
}

\author{Yingying Wei, Hongkai Ji}

\references{Yingying Wei, Xia Li, Qianfei Wang, Hongkai Ji(2012) iASeq: integrating multiple ChIP-seq datasets for detecting allele-specific binding.}

\examples{
data(sampleASE)
#fit 1 to 2 non-null correlation motifs to the data
motif.fitted<-iASeqmotif(sampleASE_exprs,sampleASE_studyid,sampleASE_repid,sampleASE_refid,
	K=1:2,iter.max=2,tol=1e-3)
}
