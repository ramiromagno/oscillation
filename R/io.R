#' Reads a GEO GSM (.csv or .csv.gz) file
#'
#' This function reads a GEO GSM file with a matrix of expression data. Rows are
#' genes (or isoforms) and columns are samples. It expects a first row of column
#' names. The first column is then used to define the rownames of the matrix.
#'
#' @param filename The path to the GSM file. It can read compressed (.gz) files.
#' 
#' @return A matrix of gene expression values. Rows are genes/isoforms and
#'   columns are samples. Rownames are gene names.
#' @export 
read_gsm_as_matrix <- function(filename) {
  dt <- data.table::fread(input = filename)
  m <- as.matrix(dt[,-1])
  rownames(m) <- dt[[1]]
  return(m)
}