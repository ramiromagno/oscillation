% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_statistics.R
\name{median2}
\alias{median2}
\title{Median over positive values}
\usage{
median2(x, na_rm = FALSE)
}
\arguments{
\item{x}{A numeric vector.}

\item{na_rm}{A logical value indicating whether NA values should be stripped
before the computation proceeds.}
}
\value{
The median value, a scalar.
}
\description{
This function computes the median over strictly positive values, \code{Inf}
excluded.
}
