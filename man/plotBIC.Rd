\name{plotBIC}
\alias{plotBIC}
\title{BIC Plot}
\description{This function plots BIC values for all fitted motif models.}

\usage{
plotBIC(fitted_cormotif)
}

\arguments{
\item{fitted_cormotif}{	The object obtained from iASeq. }
}

\author{Yingying Wei}

\seealso{
\code{\link{iASeqmotif}}, \code{\link{plotMotif}}, \code{\link{sampleASE}}
}

\examples{
example(iASeqmotif) # compute 'motif.fitted'
plotBIC(motif.fitted)
}

