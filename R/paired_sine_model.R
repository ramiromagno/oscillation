#' Paired-sine model squared error term
#' 
#' This function returns the squared error term as described in
#' \href{https://doi.org/10.1038/nmeth.3549}{Leng et al. 2015}, section
#' \emph{Oscope: paired-sine model}.
#' 
#' The definition of the error term is:
#' \deqn{\epsilon_{1,2}^2 = \sum_s ( x_{1,s}^2 + x_{2,s}^2 - 2 x_{1,s} x_{2,s} cos(\psi) - sin(\psi)^2)^2}
#' 
#' @param psi Phase shift between the gene expression of the two genes.
#' @param x1 Vector of gene expression values for gene 1.
#' @param x2 Vector of gene expression values for gene 2.
#' 
#' @return Squared error term, i.e., \eqn{\epsilon_{1,2}^2}.
#' 
#' @export
squared_epsilon <- Vectorize(
  function(psi, x1, x2) {
    sum((x1 ^ 2 + x2 ^ 2 - 2 * x1 * x2 * cos(psi) - sin(psi) ^ 2) ^ 2)
    },
  'psi'
  )

#' Optimal phase shift
#' 
#' This function determines the optimal phase shift between two gene expression
#' profiles \code{x1} and \code{x2} according to the paired-sine model, see
#' \href{https://doi.org/10.1038/nmeth.3549}{Leng et al. 2015}, section
#' \emph{Oscope: paired-sine model}.
#' 
#' @param x1 Vector of gene expression values for gene 1.
#' @param x2 Vector of gene expression values for gene 2.
#' @param psi_start The phase shift (psi) starting value for the optimisation
#'   search.
#' @param method The name of the optimisation method to be used by
#'   \code{\link[stats]{optim}}.
#' 
#' @return A \code{\link[tibble]{tibble}} where rows are for each pair of genes,
#'   and five columns: gene 1, gene 2, \code{psi}, \code{se_epsilon} and
#'   \code{score}.
#' 
#' @export
find_optimal_phase_shift <- function(x1, x2, psi_start = 0., method = "L-BFGS-B") {
  
  out <- optim(
    par = psi_start,
    fn = squared_epsilon,
    method = method,
    lower = 0.,
    upper = 2. * pi,
    x1 = x1,
    x2 = x2
  )
  
  tibble::tibble(psi = out$par,
                 sq_epsilon = out$value,
                 score = -log10(out$value)
                 )
}

#' Perform paired-sine mode analysis
#'
#' Following the paired-sine model, as described in
#' \href{https://doi.org/10.1038/nmeth.3549}{Leng et al. 2015}, section
#' \emph{Oscope: paired-sine model}, this function estimates the optimal phase
#' shift for all gene pairs by minimising the squared epsilon error term, i.e.,
#' \eqn{\epsilon_{1,2}^2}.
#'
#' @param m Matrix of gene expression values. Rows are genes and columns are samples.
#' @param verbose Whether to be chatty.
#' @param sort_by_score Whether to descending sort the final dataframe by the
#'   score column, i.e., by \eqn{-log10(sq_epsilon)}. If \code{FALSE}, the
#'   returned dataframe will be sort by the original gene pair combinations.
#' @param parallel Whether to split the computation across the virtual cores.
#' @param cores The number of virtual cores (vCPUs) to use if \code{parallel}
#'   equals \code{TRUE}.
#'
#' @return A dataframe with the optimal phase estimates. Each row corresponds to
#'   a gene pair combination. Columns:
#' \describe{
#'   \item{gene_1}{Name of gene 1.}
#'   \item{gene_2}{Name of gene 2.}
#'   \item{psi}{Optimal phase shift.}
#'   \item{sq_epsilon}{Squared epsilon (\eqn{\epsilon_{1,2}^2}).}
#'   \item{score}{Score defined as \eqn{-log10(\epsilon_{1,2}^2}). This metric
#'   can be used to sort the table of results, and select candidate oscillatory
#'   genes as those in the top gene pairs.}
#' }
#' 
#' @export
paired_sine_analysis <- function(m,
                                 verbose = FALSE,
                                 sort_by_score = FALSE,
                                 parallel = FALSE,
                                 cores = parallel::detectCores() - 2L
                                 ) {
  
  if (!is.matrix(m))
    stop("m must be a matrix.")
  
  n_row <- nrow(m)
  if (n_row < 2)
    stop("At least two genes (rows) are needed for pairwise combinations.")

  if (is.null(rownames(m)))
    stop('m must have rownames (gene names)')
  
  if (verbose) {
    message('genes: ', n_row)
    message('samples: ', ncol(m))
    message('sort_by_score: ', sort_by_score)
    message('parallel: ', parallel)
    message('number of cores: ', cores)
  }
  
  pairs_lst <- combn(seq_len(n_row), m = 2, simplify = FALSE)
  pairs_matrix <- combn(rownames(m), m = 2, simplify = TRUE)
  fn <- function(pair) find_optimal_phase_shift(m[pair[1], ], m[pair[2], ])
  
  optim_df <-
    if (parallel)
      dplyr::bind_rows(parallel::mclapply(pairs_lst, fn, mc.cores = cores))
  else
    purrr::map_dfr(.x = pairs_lst, .f = fn)

  # Add gene names  
  df <- dplyr::bind_cols(gene_1 = pairs_matrix[1,],
                         gene_2 = pairs_matrix[2,],
                         optim_df)
  
  # Sort by score?
  df <-
    if (sort_by_score)
      dplyr::arrange(df, dplyr::desc(score))
  else
    df
  
  return(df)
}
