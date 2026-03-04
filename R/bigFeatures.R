# =========================================================================
# CREATE bigPheno S3 CLASS
# =========================================================================

#' Create bigPheno object
#'
#' A list-based S3 class for storing phenotype data (e.g. gene expression)
#' in FBM format with samples on rows and phenotypes on columns.
#'
#' @param matrix Numeric matrix with samples on rows and phenotypes on columns
#' @param rowData data.frame with sample metadata. If NULL, uses rownames to create sample_name column
#' @param colData data.frame with phenotype metadata. If NULL, uses colnames to create pheno_name column
#'
#' @return A bigPheno object (list with elements: pheno, rowData, colData)
#' @export
bigPheno <- function(matrix, rowData = NULL, colData = NULL) {
  
  # Validate input matrix
  if (!is.matrix(matrix)) {
    stop("matrix must be a numeric matrix")
  }
  
  n_samples <- nrow(matrix)
  n_phenos <- ncol(matrix)
  
  # Convert matrix to FBM
  fbm <- bigstatsr::as_FBM(matrix)
  
  # Create rowData (sample metadata) if NULL
  if (is.null(rowData)) {
    sample_names <- rownames(matrix)
    if (is.null(sample_names)) {
      sample_names <- paste0("sample_", 1:n_samples)
    }
    rowData <- data.frame(sample_name = sample_names, 
                          stringsAsFactors = FALSE)
  } else {
    if (!inherits(rowData, "data.frame")) {
      stop("rowData must be a data.frame")
    }
    if (nrow(rowData) != n_samples) {
      stop(sprintf("rowData has %d rows but matrix has %d rows",
                   nrow(rowData), n_samples))
    }
  }
  
  # Create colData (phenotype metadata) if NULL
  if (is.null(colData)) {
    pheno_names <- colnames(matrix)
    if (is.null(pheno_names)) {
      pheno_names <- paste0("pheno_", 1:n_phenos)
    }
    colData <- data.frame(pheno_name = pheno_names,
                          stringsAsFactors = FALSE)
  } else {
    if (!inherits(colData, "data.frame")) {
      stop("colData must be a data.frame")
    }
    if (nrow(colData) != n_phenos) {
      stop(sprintf("colData has %d rows but matrix has %d columns",
                   nrow(colData), n_phenos))
    }
    if (!("pheno_name" %in% colnames(colData))) {
      stop("colData must contain a 'pheno_name' column")
    }
  }
  
  # Create object
  obj <- list(
    pheno = fbm,
    rowData = rowData,
    colData = colData
  )
  
  class(obj) <- "bigPheno"
  return(obj)
}
