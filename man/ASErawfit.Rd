\name{ASErawfit}
\alias{ASErawfit}
\title{Single Study Based Statistics for Allele Specific Events}
\description{This function produces standard statistics for allele-specific events based on a single RNAseq or ChIPseq study. It first pools replicates within a given study to sum the read counts for the reference allele and the non-reference allele. Then based on the pooled read counts, it calculates naive z statistic, naive Bayes statistic and empirical Bayes statistic. }

\usage{
ASErawfit(exprs,studyid,repid,refid)
}

\arguments{
\item{exprs}{ A matrix, each row of the matrix corresponds to a heterozygotic SNP and each column of the matrix corresponds to the reads count for either the reference allele or non-reference allele in a replicate of a study. }
\item{studyid}{ The group label for each column of exprs matrix. all columns in the same study have the same studyid. }
\item{repid}{ The sample label for each column of exprs matrix. The two columns within the same sample, one for reference allele and the other for non-reference allele, have the same repid. In other words, repid discriminates the different replicates within the same study.}
\item{refid}{ The reference allele label for each column of exprs matrix. Please code 0 for reference allele columns and 1 for non-reference allele columns to make the interpretation of over expressed (or bound) to be skewing to the reference allele. Otherwise, just interpret the other way round.}
}

\details{One should indicate the studyid, repid and refid for each column clearly.
}

\value{
\item{z}{ Naive z statistic. A matrix, each row of the matrix corresponds to a heteroygpotic SNP of the input matrix ('exprs') and each column corresponds to a study.}
\item{b}{ Naive Bayes statistic. A matrix, each row of the matrix corresponds to a heteroygpotic SNP of the input matrix ('exprs') and each column corresponds to a study.}
\item{B}{ Empirical Bayes statistic. A matrix, each row of the matrix corresponds to a heteroygpotic SNP of the input matrix ('exprs') and each column corresponds to a study.}
\item{c0d}{ \eqn{\alpha} parameter for the null beta prior distribution for pooled counts for each study. A vector whose length equals to the number of studies.}
\item{d0d}{ \eqn{\beta} parameter for the null beta prior distribution for pooled counts for each study. A vector whose length equals to the number of studies.}
\item{p0d}{ Mean of the null beta prior distribution for pooled counts for each study. A vector whose length equals to the number of studies.}
\item{p0dz}{ Raw mean of the reference allele proportion. A vector whose length equals to the number of studies.}

}

\author{Yingying Wei}


\references{Yingying Wei, Xia Li, Qianfei Wang, Hongkai Ji (2012) iASeq: integrating multiple ChIP-seq datasets for detecting allele-specific binding.}

\seealso{
\code{\link{sampleASE}} 
}
\examples{
data(sampleASE)
raw.fitted<-ASErawfit(sampleASE_exprs,sampleASE_studyid,sampleASE_repid,sampleASE_refid)
}
