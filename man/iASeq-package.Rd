\name{iASeq-package}
\alias{iASeq-package}
\alias{iASeq}
\docType{package}
\title{
iASeq: integrating multiple sequencing datasets for detecting allele-specific events
}
\description{
In diploid organisms, certain genes can be expressed, methylated or regulated in an allele-specific manner, corresponding to allele-specific expression, allele-specific methylation and allele-specific binding. These allele-specific events (AS) are of high interest for phenotypic diversity and disease susceptibility. Next generation sequencing technologies provide opportunities to study AS globally. However, little is known about the mechanism of AS. For instance, the patterns of allele-specific binding across different Transcription Factors (TFs) and histone modifications (HMs) are unclear. Moreover, the limited number of reads on heterozygotic SNPs results in low-signal-to-noise ratio when calling AS. Here, we propose a Bayes hierarchical model to study AS by jointly analyzing multiple ChIPseq studies, RNAseq studies or MeDIPseq studies. The model is able to learn the patterns of AS across studies and make substantial improvement in calling AS.}
\details{
\tabular{ll}{
Package: \tab iASeq\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2012-02-13\cr
License: \tab GPL-2\cr
}

}
\author{
Yingying Wei, Hongkai Ji

Maintainer: Yingying Wei <ywei@jhsph.edu>
}
\references{Yingying Wei, Xia Li, Qianfei Wang, Hongkai Ji (2012) iASeq: integrating multiple ChIP-seq datasets for detecting allele-specific binding.}

\seealso{
\code{\link{iASeqmotif}}, \code{\link{plotBIC}}, \code{\link{plotMotif}}, \code{\link{sampleASE}}, \code{\link{ASErawfit}}, \code{\link{singleEMfit}}, \code{\link{sampleASE}}
}

