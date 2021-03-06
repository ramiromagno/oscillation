% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score_matrix.R
\name{score_matrix}
\alias{score_matrix}
\title{Converts the data frame of paired-sine scores to a matrix}
\usage{
score_matrix(df, gene_1 = "gene_1", gene_2 = "gene_2", score = "score")
}
\arguments{
\item{df}{A data frame as returned by
\code{\link[oscillation]{paired_sine_analysis}}.}

\item{gene_1}{Name of the column containing the names of the first gene in
the gene pair (default is \code{"gene_1"}).}

\item{gene_2}{Name of the column containing the names of the second gene in
the gene pair (default is \code{"gene_2"}).}

\item{score}{Name of the column containing the scores (default is
\code{"score"}).}
}
\value{
A symmetric matrix of gene pair scores
  (\eqn{-log10(\epsilon_{1,2}^2})). Column and row names are the gene names
  in alphabetical order.
}
\description{
This function expects a data frame \code{df}, as returned by the function
\code{\link[oscillation]{paired_sine_analysis}}, i.e., containing the scores
(\eqn{-log10(\epsilon_{1,2}^2})) for each gene pair. The resulting matrix is
a wide format of the gene-pairs scores (originally in long format in
\code{df}).
}
\examples{
m <- cbind(t(combn(1:5, m = 2, simplify = TRUE)), 1:10)
colnames(m) <- c('x', 'y', 'z')

score_matrix(m, gene_1 = 'x', gene_2 = 'y', score = 'z')
 
}
