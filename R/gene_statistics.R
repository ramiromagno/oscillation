#' Gene expression statistics
#' 
#' This function computes five statistics: mean, median, variance, \eqn{q}
#' (actually, \code{q1} and \code{q2}), and \eqn{\phi}. The statistics \eqn{q}
#' and \eqn{\phi} are parameter estimations related to the
#' \href{https://en.wikipedia.org/wiki/Negative_binomial_distribution}{Negative
#' Binomial} distribution, a commonly assumed model for the distribution of gene
#' (or isoform) expression counts.
#' 
#' @param m Matrix of gene expression values. Rows are genes and columns are
#'   samples.
#' @param p_min Minimum value for parameter p of the Negative Binomial
#'   distribution to be assumed when estimating q.
#' 
#' @export
gene_statistics <- function(m, p_min = 0) {
  
  if (!is.matrix(m))
    stop("m must be a matrix.")
  
  gene_names <- rownames(m)
  if (is.null(gene_names))
    stop('m must have rownames (gene names)')
  
  means <- rowMeans(m)
  medians <- apply(m, 1, median)
  variances <- apply(m, 1, var)
  
  # Negative binomial parameter q: q = 1 - p
  q1 <- means / variances
  
  # Here we cap the estimation of q to a value close to 1.
  q_max <- 1 - p_min
  
  # q2 is q1 after applying the cap q_max.
  q2 <- pmin(q1, q_max)
  
  # The estimation of phi is here based on q2
  # instead of q1.
  phi <- (1 - q2) / (means * q2)
  
  df <- tibble::tibble(
    gene = gene_names,
    mean = means,
    median = medians,
    variance = variances,
    q1 = q1,
    q2 = q2,
    phi = phi
  )
  
  return(df)
}

#' Filter genes by mean counts
#' 
#' This functions takes a dataframe of gene expression statistics (see function
#' \code{\link[oscillation]{gene_statistics}}) and keeps only those genes whose
#' mean is within \code{low} and \code{high}, i.e., lie in the interval
#' [\code{low}, \code{high}].
#' 
#' @param gene_stat A dataframe as returned by function
#'   \code{\link[oscillation]{gene_statistics}}).
#' @param low Lower bound for the mean (use \code{-Inf} for an effective no
#'   bound).
#' @param high Upper bound for the mean (use \code{Inf} for an effective no
#'   bound).
#' 
#' @return A dataframe of the same structure as \code{gene_stat} of filtered
#'   genes.
#' 
#' @export
filter_by_mean_count <- function(gene_stat, low = 100, high = Inf) {
  
  gene_stat %>%
    dplyr::filter(mean >= low & mean <= high)
  
}

#' Linear regression on gene variance and mean
#'
#' This function performs a linear regression on the log-transformed values of
#' the variance vs mean, i.e., \code{log(variance)~log(mean)}. Genes whose
#' (log(mean), log(variance)) points stay above the fitted line can be
#' considered high variance genes. These values will have a \strong{negative}
#' residual (\code{.resid} column in dataframe \code{residuals}).
#'
#' To keep only the high variance genes, you can filter by residual:
#' \code{mean_variance_fit(foo)$residuals \%>\% dplyr::filter(.resid < 0)}.
#' 
#' @param gene_stat_df A dataframe of at least three columns: \code{gene},
#'   \code{mean}, and \code{variance}.
#'
#' @return A list of three dataframes:
#' \describe{
#'   \item{regression_coeficients}{Information about the regression coeficients. Each row is a regression coefficient.}
#'   \item{goodness_of_fit}{Information on the goodness of fit.}
#'   \item{residuals}{Information about each fitted point. The column
#'   \code{.resid} contains the difference between the fitted line and the data
#'   point, i.e., \eqn{\hat{y}-y}. Please note this is not the most customary
#'   definition of residual:
#'   \href{https://github.com/tidymodels/broom/issues/802}{broom issue 802}.
#'   But, given this definition, the high variance genes are indeed those that
#'   show a \strong{negative} residual.}
#' }  
#'
#' @export
mean_variance_fit <- function(gene_stat_df) {
  
  df <- gene_stat_df[c('gene', 'mean', 'variance')]
  
  fit <- lm(log(variance) ~ log(mean), data = df)
  
  list(
    regression_coeficients = broom::tidy(fit),
    goodness_of_fit = broom::glance(fit),
    residuals = broom::augment(fit, data = df)
    )
  
}

#' Calculate the pseudo-reference sample
#' 
#' This function calculates the geometric mean across samples. This can be
#' regarded as a pseudo-reference sample. See
#' \href{http://genomebiology.com/2010/11/10/R106}{Anders & Huber (2010)}, page
#' 3, denominator in equation 5.
#' 
#' @param m Matrix of gene expression values. Rows are genes and columns are
#'   samples.
#'   
#' @return A vector of geometric means, one for each gene (can be regarded as
#'   pseudo-reference sample for each gene).
#'   
#' @export   
pseudo_ref_sample <- function(m) exp(rowMeans(log(m)))

#' Median normalisation size factors
#' 
#' This function calculates the median normalisation size factors of a matrix of gene
#' expression counts according to
#' \href{http://genomebiology.com/2010/11/10/R106}{Anders & Huber (2010)}, page
#' 3, equation 5.
#' 
#' @param m Matrix of gene expression values. Rows are genes and columns are
#'   samples.
#'
#' @return A vector of size factors, one for each sample.
#'
#' @export
median_norm_size_factors <- function(m) {
  
  geom_means <- pseudo_ref_sample(m)
  # matrixStats::colMedians((m / geom_means)[geom_means > 0, , drop = FALSE])
  apply(m, 2, function(j) median((j/geom_means)[geom_means > 0]))
}

#' Median over positive values
#' 
#' This function computes the median over strictly positive values, \code{Inf}
#' excluded.
#' 
#' @param x A numeric vector.
#'
#' @return The median value, a scalar.
#'
#' @export
median2 <- function(x) median(x[x > 0 & x < Inf])

#' Geometric mean
#' 
#' This function computes the geometric mean.
#' 
#' @param x A numeric vector.
#' @param na.rm Whether to remove \code{NA} before calculation.
#' @param zero.propagate Whether zeros are allowed to propagate in the
#'   calculation or whether they are treated as \code{NA}.
#'
#' @return The geometric mean.
#'
#' @export
gm_mean = function(x,
                   na_rm = TRUE,
                   zero_propagate = TRUE) {
  if (any(x < 0, na.rm = TRUE)) {
    return(NaN)
  }
  if (zero_propagate) {
    if (any(x == 0, na.rm = TRUE)) {
      return(0)
    }
    exp(mean(log(x), na.rm = na_rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm = na_rm) / length(x))
  }
}

#' Median normalisation size factors (alternative)
#' 
#' This function calculates the median normalisation size factors of a matrix of
#' gene expression counts following the code of the function
#' \code{\link[EBSeq]{MedianNorm}} when the parameter \code{alternative=TRUE}.
#' 
#' @param m Matrix of gene expression values. Rows are genes and columns are
#'   samples.
#'
#' @return A vector of size factors, one for each sample.
#'
#' @export
median_norm_size_factors2 <- function(m) {

  # m2: a square matrix whose dimension is the same as the number of columns of
  # m.
  m2 <- apply(m, 2, function(j) apply(m, 2, function(jj) median2(j/jj)))

  sample_indices <- seq_len(ncol(m))
  
  # The index combination -j, j ensures that the diagonal values of m2 are
  # excluded from the computation.
  sapply(sample_indices, function(j) gm_mean(m2[-j, j]))
  
}