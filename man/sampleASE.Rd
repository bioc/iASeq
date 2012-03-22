\name{sampleASE}
\alias{sampleASE}
\alias{sampleASE_exprs}
\alias{sampleASE_refid}
\alias{sampleASE_repid}
\alias{sampleASE_studyid}
\title{Example Dataset for iASeq}
\description{Here we present four files needed for the various iASeq fit functions.}

\details{sampleASE consists of five ChIP-seq studies from ENCODE GM12878 cell lines with 5504 heterozygotic SNPs. Each study has two replicates. Each replicate's fastq reads file was aligned to hg18 whole genome using MAQ (Version 0.7.1) with default parameters. Uniquely alignments were extracted following the mapping quality above 0. Alignment can also be down using other alignment tools such as Bowtie. The GM12878 genotype data was downloaded from the website http://alleleseq.gersteinlab.org/downloads.html [Rozowsky J et al.]. The reads aligned to each allele of a heterozygotic SNP were counted correspondingly.  sampleASE_exprs saves the read counts. sampleASE_studyid prepares the study label for each sample; sample_repid describes the sample label for each column; sample_refid shows whether each column corresponds to the reference allele or the non-reference allele.}

\value{
\item{sampleASE_exprs}{ The read count matrix for the example dataset used by iASeq package. Each row of the matrix corresponds to a heterozygotic SNP and each column of the matrix corresponds to the reads count for either the reference allele or non-reference allele in a replicate of a study. }
\item{sampleASE_studyid}{ The group label for each column of sampleASE_exprs matrix. All columns in the same study have the same studyid and there are five ChIP-seq studies in this example. }
\item{sampleASE_repid}{ The sample label for each column of sampleASE_exprs matrix. The two columns within the same sample, one for reference allele and the other for non-reference allele, have the same repid. In other words, repid discriminates the different replicates within the same study. Here each study has two replicates.}
\item{sampleASE_refid}{ The reference allele label for each column of sampleASE_exprs matrix. 0 is coded for reference allele columns and 1 is coded for non-reference allele columns.}
}

\references{Yingying Wei, Xia Li, Qianfei Wang, Hongkai Ji (2012) iASeq: integrating multiple ChIP-seq datasets for detecting allele-specific binding.
Rozowsky J, Abyzov A, Wang J, Alves P, Raha D, Harmanci A, Leng J, Bjornson R, Kong Y, Kitabayashi N,
et al (2011) AlleleSeq: analysis of allele-specific expression and binding in a network framework. Mol Syst
Biol 7:522.
}

\seealso{
\code{\link{iASeqmotif}}, \code{\link{plotBIC}}, \code{\link{plotMotif}}, \code{\link{ASErawfit}}, \code{\link{singleEMfit}}
}
\examples{
data(sampleASE)
#fit 1 to 2 non-null correlation motifs to the data
motif.fitted<-iASeqmotif(sampleASE_exprs,sampleASE_studyid,sampleASE_repid,sampleASE_refid,
	K=1:2,iter.max=2,tol=1e-3)
plotBIC(motif.fitted)
plotMotif(motif.fitted$bestmotif,cutoff=0.9)
}