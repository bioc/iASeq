\name{singleEMfit}
\alias{singleEMfit}
\title{Single Study EM Fit for Allele Specific Events}
\description{This function runs an EM algorithm for allele-specific events based on a single RNAseq or ChIPseq study. It first pools replicates within a given study to sum the read counts for the reference allele and the non-reference allele. Then based on the pooled read counts, it fits an EM algorithm with three mixture components, the null distribution, the reference allele over expressed (bound) and under expressed (bound) distributions to the data.}

\usage{
singleEMfit(exprs,studyid,repid,refid,iter.max=100,tol=1e-3)
}

\arguments{
\item{exprs}{ A matrix, each row of the matrix corresponds to a heterozygotic SNP and each column of the matrix corresponds to the reads count for either the reference allele or non-reference allele in a replicate of a study. }
\item{studyid}{ The group label for each column of exprs matrix. All columns in the same study have the same studyid. }
\item{repid}{ The sample label for each column of exprs matrix. The two columns within the same sample, one for reference allele and the other for non-reference allele, have the same repid. In other words, repid discriminates the different replicates within the same study.}
\item{refid}{ The reference allele label for each column of exprs matrix. Please code 0 for reference allele columns and 1 for non-reference allele columns to make the interpretation of over expressed(or bound) to be skewing to the reference allele. Otherwise, just interpret the other way round.}
\item{tol}{ The relative tolerance level of error. }
\item{iter.max}{ Maximun number of iterations. }
}

\value{
\item{p.study}{ The posterior probability for each SNP to be allele-specific event within each study. A matrix where each row corresponds to a SNP and each column corresponds to a study.}
\item{motif.qup}{ Fitted values of probability for the reference allele of each SNP to be over expressed (or bound) within each study. A matrix where each row corresponds to a SNP and each column corresponds to a study.}
\item{motif.qdown}{ Fitted values of probability for  the reference allele of each SNP to be under expressed (or bound) within each study. A matrix where each row corresponds to a SNP and each column corresponds to a study.}
\item{condlike}{  A list where each element is a matrix and corresponds to a study. Each row of each matrix corresponds to a SNP. The three column of each matrix represents the posterior probability for a SNP to belong to the null distribution, the over expressed distribution and the under expressed distribution within the given study.}
}

\author{Yingying Wei}

\references{Yingying Wei, Xia Li, Qianfei Wang, Hongkai Ji (2012) iASeq: integrating multiple ChIP-seq datasets for detecting allele-specific binding.}

\seealso{
\code{\link{sampleASE}}
}
\examples{
data(sampleASE)
singleEM.fitted<-singleEMfit(sampleASE_exprs,sampleASE_studyid,sampleASE_repid,
		sampleASE_refid,iter.max=2,tol=1e-3)
}
