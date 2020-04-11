#' Converts the data frame of paired-sine scores to a matrix
#'
#' This function expects a data frame \code{df}, as returned by the function
#' \code{\link[oscillation]{paired_sine_analysis}}, i.e., containing the scores
#' (\eqn{-log10(\epsilon_{1,2}^2})) for each gene pair. The resulting matrix is
#' a wide format of the gene-pairs scores (originally in long format in
#' \code{df}).
#'
#' @param df A data frame as returned by
#'   \code{\link[oscillation]{paired_sine_analysis}}.
#' @param gene_1 Name of the column containing the names of the first gene in
#'   the gene pair (default is \code{"gene_1"}).
#' @param gene_2 Name of the column containing the names of the second gene in
#'   the gene pair (default is \code{"gene_2"}).
#' @param score Name of the column containing the scores (default is
#'   \code{"score"}).
#' 
#' @return A symmetric matrix of gene pair scores
#'   (\eqn{-log10(\epsilon_{1,2}^2})). Column and row names are the gene names
#'   in alphabetical order.
#'   
#' @examples
#' m <- cbind(t(combn(1:5, m = 2, simplify = TRUE)), 1:10)
#' colnames(m) <- c('x', 'y', 'z')
#' 
#' score_matrix(m, gene_1 = 'x', gene_2 = 'y', score = 'z')
#'  
#' @export 
score_matrix <- function(df,
                         gene_1 = 'gene_1',
                         gene_2 = 'gene_2',
                         score = 'score') {
  
  # genes: vector of all the gene names.
  genes <- sort(union(df[, gene_1], df[, gene_2]))
  
  # n_genes: number of unique gene names.
  n_genes <- length(genes)
  
  # m: Allocating a symmetric matrix of genes by genes.
  m <-
    matrix(
      data = NA_real_,
      nrow = n_genes,
      ncol = n_genes,
      dimnames = list(genes, genes)
    )
  
  # Populate the matrix m with the scores, and use subsetting by name
  # to indicate which score goes to which (i, j) coordinate in m.
  m[cbind(df[, gene_1], df[, gene_2])] <- df[, score]
  #m[cbind(df[, gene_2], df[, gene_1])] <- df[, score]
  
  # Change NAs in lower triangular matrix to values found in upper triangular
  # matrix
  na.index <- which(is.na(m), arr.ind = TRUE)
  m[na.index] <- m[na.index[, c(2, 1)]]
  
  # Set the diagonal values to zero.
  m[cbind(genes, genes)] <- 0.
  
  return(m)
}
