% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PQPQ_Mat.R
\name{pqpq_orig}
\alias{pqpq_orig}
\title{Matlab's PQPQ function}
\usage{
pqpq_orig(pepDataFixed, confList)
}
\arguments{
\item{pepDataFixed}{An object of type 'pepObject' (from pmartRqc), which has been preferably
normalized and filtered. The matrix 'e_data' has n rows.}

\item{confList}{A list of n confidence values from 0 to 100, corresponding to the peptides in 'e_data'.}
}
\value{
A data frame with 3 columns: \describe{
\item{Protein}{Which protein this peptide is associated with}
\item{Peptide}{The row name of the valid peptide}
\item{ID}{Which cluster the peptide is considered clustered to}}
}
\description{
This function implements the Matlab version of PQPQ, removing all peptides with missing values and
performing PQPQ on the remaining data.
}
\seealso{
\code{\link{pqpq}}
}
\keyword{Matlab}
\keyword{PQPQ}
