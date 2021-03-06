% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_statistics.R
\name{median_norm_size_factors}
\alias{median_norm_size_factors}
\title{Median normalisation size factors}
\usage{
median_norm_size_factors(m)
}
\arguments{
\item{m}{Matrix of gene expression values. Rows are genes and columns are
samples.}
}
\value{
A vector of size factors, one for each sample.
}
\description{
This function calculates the median normalisation size factors of a matrix of gene
expression counts according to
\href{http://genomebiology.com/2010/11/10/R106}{Anders & Huber (2010)}, page
3, equation 5.
}
