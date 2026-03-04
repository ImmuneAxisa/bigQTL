# =========================================================================
# MAIN FUNCTION: Conditional QTL analysis
# =========================================================================

#' Run conditional QTL analysis
#'
#' Performs conditional cis-QTL analysis using stepwise conditioning and
#' all-but-one conditioning on file-backed genotype and phenotype matrices.
#'
#' @param bigpheno bigPheno object with phenotype data (samples x phenotypes)
#' @param bigsnp bigSNP object with genotype data (samples x snps)
#' @param pheno_coord Data frame with pheno_name, chromosome, start, end columns
#' @param design_base Data frame with covariates; row names = sample IDs
#' @param cis_window Padding around phenotype coordinates (default 1e6)
#' @param do_conditioning Logical; perform stepwise conditioning (default TRUE)
#' @param pval_threshold P-value threshold for stepwise conditioning (default 1e-3)
#' @param do_allbutone Logical; perform all-but-one conditioning (default TRUE)
#' @param do_rint Logical; apply RINT transformation to each phenotype (default TRUE)
#' @param ncores Cores for per-SNP regression within-phenotype (default 1)
#' @param ncores_phenos Cores for across-phenotype parallelisation (default 1)
#' @param output_dir Directory for Parquet output (default "./qtl_results")
#'
#' @return List with elements: stepwise (Arrow dataset), allbutone (Arrow dataset or NULL),
#'   stepwise_dir, allbutone_dir
#' @export
run_conditional_qtl <- function(
    bigpheno,                  # bigPheno object with phenotype data (samples x phenotypes)
    bigsnp,                    # bigSNP object with genotype data (samples x snps)
    pheno_coord,               # Data frame with pheno_name, chromosome, start, end columns
    design_base,               # Data frame with covariates; row names = sample IDs
    cis_window = 1e6,          # Padding around phenotype coordinates
    do_conditioning = TRUE,
    pval_threshold = 1e-3,
    do_allbutone = TRUE,
    do_rint = TRUE,            # Apply RINT transformation to each phenotype (default TRUE)
    ncores = 1,                # Cores for per-SNP regression (within-phenotype)
    ncores_phenos = 1,         # Cores for across-phenotype parallelisation
    output_dir = "./qtl_results") {
  
  # =====================================================================
  # Resolve sample IDs from design_base row names
  # =====================================================================
  
  sample_ids <- rownames(design_base)
  
  if (is.null(sample_ids) || identical(sample_ids, as.character(seq_len(nrow(design_base))))) {
    stop("design_base must have meaningful row names (sample IDs), not default integer row names")
  }
  
  # Map to genotype row indices
  ind.row.snp <- match(sample_ids, bigsnp$fam$sample.ID)
  if (any(is.na(ind.row.snp))) {
    missing <- sample_ids[is.na(ind.row.snp)]
    stop(sprintf("%d sample ID(s) from design_base not found in bigsnp$fam$sample.ID: %s",
                 length(missing), paste(head(missing, 3), collapse = ", ")))
  }
  
  # Map to phenotype row indices
  ind.row.pheno <- match(sample_ids, bigpheno$rowData$sample_name)
  if (any(is.na(ind.row.pheno))) {
    missing <- sample_ids[is.na(ind.row.pheno)]
    stop(sprintf("%d sample ID(s) from design_base not found in bigpheno$rowData$sample_name: %s",
                 length(missing), paste(head(missing, 3), collapse = ", ")))
  }
  
  message(sprintf("  Matched %d samples from design_base to bigsnp and bigpheno",
                  length(sample_ids)))
  
  # Validate pheno_coord
  required_cols <- c("pheno_name", "chromosome", "start", "end")
  if (!all(required_cols %in% colnames(pheno_coord))) {
    stop(sprintf("pheno_coord must contain columns: %s",
                 paste(required_cols, collapse = ", ")))
  }
  
  # Create output directories
  stepwise_dir <- file.path(output_dir, "stepwise")
  allbutone_dir <- file.path(output_dir, "allbutone")
  
  if (!dir.exists(stepwise_dir)) dir.create(stepwise_dir, recursive = TRUE)
  if (!dir.exists(allbutone_dir)) dir.create(allbutone_dir, recursive = TRUE)
  
  # =====================================================================
  # Process phenotypes (parallelised with mclapply or sequential with lapply)
  # =====================================================================
  
  phenos <- pheno_coord$pheno_name
  
  process_fn <- function(pheno) {
    process_pheno(
      pheno          = pheno,
      bigpheno       = bigpheno,
      bigsnp         = bigsnp,
      pheno_coord    = pheno_coord,
      design_base    = design_base,
      ind.row.snp    = ind.row.snp,
      ind.row.pheno  = ind.row.pheno,
      cis_window     = cis_window,
      do_conditioning = do_conditioning,
      pval_threshold = pval_threshold,
      do_allbutone   = do_allbutone,
      do_rint        = do_rint,
      ncores         = ncores,
      stepwise_dir   = stepwise_dir,
      allbutone_dir  = allbutone_dir
    )
  }
  
  if (ncores_phenos > 1) {
    parallel::mclapply(phenos, process_fn, mc.cores = ncores_phenos)
  } else {
    lapply(phenos, process_fn)
  }
  
  # =====================================================================
  # Open and return Parquet datasets
  # =====================================================================
  
  message("\nOpening Parquet datasets...")
  
  stepwise_dataset <- arrow::open_dataset(stepwise_dir)
  
  allbutone_dataset <- if (dir.exists(allbutone_dir) && length(list.files(allbutone_dir)) > 0) {
    arrow::open_dataset(allbutone_dir)
  } else {
    NULL
  }
  
  message(sprintf("Results saved to: %s", output_dir))
  
  return(list(
    stepwise     = stepwise_dataset,
    allbutone    = allbutone_dataset,
    stepwise_dir = stepwise_dir,
    allbutone_dir = allbutone_dir
  ))
}


# =========================================================================
# WRAPPER FUNCTION: bigQTL — standard preparation + QTL analysis
# =========================================================================

#' Run bigQTL: standard preparation and conditional QTL analysis
#'
#' A convenience wrapper around \code{run_conditional_qtl()} that first
#' computes phenotype PCs and genotype PCs using the package helper
#' functions, appends them to \code{design_base}, and then runs the full
#' conditional QTL analysis.
#'
#' @param bigpheno bigPheno object with phenotype data (samples x phenotypes)
#' @param bigsnp bigSNP object with genotype data (samples x snps)
#' @param pheno_coord Data frame with pheno_name, chromosome, start, end columns
#' @param design_base Data frame with covariates; row names = sample IDs
#' @param n_pheno_pcs Number of phenotype PCs to compute and add as covariates
#'   (default 5)
#' @param n_geno_pcs Number of genotype PCs to compute and add as covariates
#'   (default 5)
#' @param exclude_pheno_names Character vector of phenotype names to exclude
#'   from phenotype PCA (e.g. sex-chromosome phenotypes). Default NULL.
#' @param ... Additional arguments passed to \code{run_conditional_qtl()}
#'
#' @return List with elements: stepwise (Arrow dataset), allbutone (Arrow dataset or NULL),
#'   stepwise_dir, allbutone_dir
#' @export
bigQTL <- function(bigpheno, bigsnp, pheno_coord, design_base,
                   n_pheno_pcs = 5, n_geno_pcs = 5,
                   exclude_pheno_names = NULL, ...) {

  keep_ids <- rownames(design_base)

  # Compute phenotype PCs and append to design
  pheno_pcs <- compute_pheno_pcs(
    bigpheno           = bigpheno,
    keep_ids           = keep_ids,
    n_pcs              = n_pheno_pcs,
    exclude_pheno_names = exclude_pheno_names
  )
  pheno_pcs_aligned <- pheno_pcs[keep_ids, , drop = FALSE]

  # Compute genotype PCs and append to design
  geno_pcs <- compute_geno_pcs(
    bigsnp   = bigsnp,
    keep_ids = keep_ids,
    n_pcs    = n_geno_pcs
  )
  geno_pcs_aligned <- geno_pcs[keep_ids, , drop = FALSE]

  design_augmented <- cbind(design_base, pheno_pcs_aligned, geno_pcs_aligned)

  run_conditional_qtl(
    bigpheno    = bigpheno,
    bigsnp      = bigsnp,
    pheno_coord = pheno_coord,
    design_base = design_augmented,
    ...
  )
}


# =========================================================================
# PER-PHENOTYPE FUNCTION: Process a single phenotype
# =========================================================================

#' Process a single phenotype for conditional QTL analysis
#'
#' @param pheno Phenotype name to process
#' @param bigpheno bigPheno object
#' @param bigsnp bigSNP object
#' @param pheno_coord Data frame with phenotype coordinates
#' @param design_base Data frame with covariates
#' @param ind.row.snp Row indices into genotype FBM for samples
#' @param ind.row.pheno Row indices into phenotype FBM for samples
#' @param cis_window Cis window size
#' @param do_conditioning Logical; perform stepwise conditioning
#' @param pval_threshold P-value threshold for stepwise conditioning
#' @param do_allbutone Logical; perform all-but-one conditioning
#' @param do_rint Logical; apply RINT to phenotype
#' @param ncores Number of cores
#' @param stepwise_dir Output directory for stepwise results
#' @param allbutone_dir Output directory for all-but-one results
#'
#' @return Invisibly NULL (results written to disk)
#' @keywords internal
process_pheno <- function(pheno, bigpheno, bigsnp, pheno_coord,
                          design_base, ind.row.snp, ind.row.pheno,
                          cis_window, do_conditioning, pval_threshold,
                          do_allbutone, do_rint, ncores,
                          stepwise_dir, allbutone_dir) {
  
  message(sprintf("Processing %s...", pheno))
  
  # Get phenotype coordinates and index
  pheno_row <- pheno_coord[pheno_coord$pheno_name == pheno, ]
  pheno_idx <- get_pheno_indices(bigpheno, pheno)
  
  # Extract phenotype for kept individuals
  y <- bigpheno$pheno[ind.row.pheno, pheno_idx]
  
  # RINT transform
  if (do_rint) {
    y <- rint(y)
  }
  
  # Get cis SNPs
  cis_result <- get_cis_snps(bigsnp, pheno_row$chromosome,
                             pheno_row$start, pheno_row$end, cis_window)
  snp_indices <- cis_result$indices
  cis_snps_pheno <- cis_result$names
  
  if (length(snp_indices) == 0) {
    message(sprintf("  No cis SNPs found for %s", pheno))
    return(invisible(NULL))
  }
  
  message(sprintf("  Found %d cis SNPs", length(snp_indices)))
  
  # ==== Step 0: Marginal testing ====
  
  results_step0 <- test_snps_with_indices(
    bigsnp      = bigsnp,
    y           = y,
    snp_indices = snp_indices,
    snp_names   = cis_snps_pheno,
    design_base = design_base,
    ind.row     = ind.row.snp,
    ncores      = ncores
  )
  results_step0$step <- 0L
  results_step0$conditioning_snps <- NA_character_
  results_step0 <- results_step0[, c("step", "conditioning_snps",
                                     setdiff(names(results_step0),
                                             c("step", "conditioning_snps")))]
  
  # ==== Stepwise conditioning ====
  
  stepwise_result <- run_stepwise(
    results_step0  = results_step0,
    bigsnp         = bigsnp,
    y              = y,
    snp_indices    = snp_indices,
    cis_snps_pheno = cis_snps_pheno,
    design_base    = design_base,
    ind.row.snp    = ind.row.snp,
    do_conditioning = do_conditioning,
    pval_threshold = pval_threshold,
    ncores         = ncores
  )
  
  stepwise_all <- stepwise_result$stepwise_tables
  conditioning_snps <- stepwise_result$conditioning_snps
  
  # ==== All-but-one conditioning ====
  
  allbutone_all <- run_allbutone(
    conditioning_snps = conditioning_snps,
    stepwise_tables   = stepwise_all,
    bigsnp            = bigsnp,
    y                 = y,
    snp_indices       = snp_indices,
    cis_snps_pheno    = cis_snps_pheno,
    design_base       = design_base,
    ind.row.snp       = ind.row.snp,
    do_allbutone      = do_allbutone,
    ncores            = ncores
  )
  
  # ==== Write results ====
  
  pheno_partition <- paste0("pheno=", pheno)
  
  # Stepwise
  sw_dir <- file.path(stepwise_dir, pheno_partition)
  dir.create(sw_dir, recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(do.call(rbind, stepwise_all), file.path(sw_dir, "part-0.parquet"))
  message(sprintf("    Wrote stepwise results to %s", pheno_partition))
  
  # All-but-one
  if (length(allbutone_all) > 0) {
    abo_dir <- file.path(allbutone_dir, pheno_partition)
    dir.create(abo_dir, recursive = TRUE, showWarnings = FALSE)
    arrow::write_parquet(do.call(rbind, allbutone_all), file.path(abo_dir, "part-0.parquet"))
    message(sprintf("    Wrote all-but-one results to %s", pheno_partition))
  }
  
  return(invisible(NULL))
}


# =========================================================================
# STEPWISE CONDITIONING
# =========================================================================

#' Run stepwise conditioning for a phenotype
#'
#' @param results_step0 Data frame of marginal association results (step 0)
#' @param bigsnp bigSNP object
#' @param y Numeric phenotype vector
#' @param snp_indices Integer vector of cis-SNP column indices
#' @param cis_snps_pheno Character vector of cis-SNP names
#' @param design_base Data frame of covariates
#' @param ind.row.snp Integer vector of row indices for samples in genotype FBM
#' @param do_conditioning Logical; perform stepwise conditioning
#' @param pval_threshold P-value threshold for stepwise conditioning
#' @param ncores Number of cores
#'
#' @return List with elements: stepwise_tables, conditioning_snps
#' @keywords internal
run_stepwise <- function(results_step0, bigsnp, y, snp_indices,
                         cis_snps_pheno, design_base, ind.row.snp,
                         do_conditioning, pval_threshold, ncores) {
  
  stepwise_tables <- list(results_step0)
  conditioning_snps <- character()
  
  if (!do_conditioning) {
    return(list(stepwise_tables = stepwise_tables,
                conditioning_snps = conditioning_snps))
  }
  
  # Check step 0 lead SNP
  lead <- results_step0[which.min(results_step0$pvalue), ]
  
  if (lead$pvalue >= pval_threshold) {
    return(list(stepwise_tables = stepwise_tables,
                conditioning_snps = conditioning_snps))
  }
  
  message(sprintf("  Step 0: Lead SNP %s passes (p = %.2e)",
                  lead$snp, lead$pvalue))
  
  conditioning_snps <- lead$snp
  step <- 1
  
  while (TRUE) {
    
    results_step <- test_snps_with_indices(
      bigsnp           = bigsnp,
      y                = y,
      snp_indices      = snp_indices,
      snp_names        = cis_snps_pheno,
      design_base      = design_base,
      ind.row          = ind.row.snp,
      snp_conditioning = conditioning_snps,
      ncores           = ncores
    )
    results_step$step <- step
    results_step$conditioning_snps <- paste(conditioning_snps, collapse = ";")
    results_step <- results_step[, c("step", "conditioning_snps",
                                     setdiff(names(results_step),
                                             c("step", "conditioning_snps")))]
    
    stepwise_tables[[length(stepwise_tables) + 1]] <- results_step
    
    new_lead <- results_step[which.min(results_step$pvalue), ]
    
    if (new_lead$pvalue < pval_threshold) {
      message(sprintf("    Step %d: New lead SNP %s passes (p = %.2e)",
                      step, new_lead$snp, new_lead$pvalue))
      conditioning_snps <- c(conditioning_snps, new_lead$snp)
      step <- step + 1
    } else {
      message(sprintf("    Step %d: Stopping (min p = %.2e)",
                      step, new_lead$pvalue))
      break
    }
  }
  
  return(list(stepwise_tables = stepwise_tables,
              conditioning_snps = conditioning_snps))
}


# =========================================================================
# ALL-BUT-ONE CONDITIONING
# =========================================================================

#' Run all-but-one conditioning for a phenotype
#'
#' @param conditioning_snps Character vector of independent SNP names
#' @param stepwise_tables List of data frames from stepwise conditioning
#' @param bigsnp bigSNP object
#' @param y Numeric phenotype vector
#' @param snp_indices Integer vector of cis-SNP column indices
#' @param cis_snps_pheno Character vector of cis-SNP names
#' @param design_base Data frame of covariates
#' @param ind.row.snp Integer vector of row indices for samples in genotype FBM
#' @param do_allbutone Logical; perform all-but-one conditioning
#' @param ncores Number of cores
#'
#' @return List of data frames with all-but-one results
#' @keywords internal
run_allbutone <- function(conditioning_snps, stepwise_tables, bigsnp, y,
                          snp_indices, cis_snps_pheno, design_base,
                          ind.row.snp, do_allbutone, ncores) {
  
  if (!do_allbutone || length(conditioning_snps) <= 1) {
    return(list())
  }
  
  message(sprintf("  Running all-but-one conditioning for %d SNPs",
                  length(conditioning_snps)))
  
  allbutone_tables <- list()
  
  for (i in seq_along(conditioning_snps)) {
    
    snps_condition_on <- conditioning_snps[-i]
    
    if (i == length(conditioning_snps)) {
      
      # Last independent SNP: reuse final stepwise step
      result <- stepwise_tables[[length(stepwise_tables)]]
      result$indep <- i
      result$conditioning_snps <- paste(snps_condition_on, collapse = ";")
      result <- result[, c("snp", "beta", "se", "t_stat", "pvalue", "fdr",
                           "indep", "conditioning_snps")]
      
    } else {
      
      result <- test_snps_with_indices(
        bigsnp           = bigsnp,
        y                = y,
        snp_indices      = snp_indices,
        snp_names        = cis_snps_pheno,
        design_base      = design_base,
        ind.row          = ind.row.snp,
        snp_conditioning = snps_condition_on,
        ncores           = ncores
      )
      result$indep <- i
      result$conditioning_snps <- paste(snps_condition_on, collapse = ";")
      result <- result[, c("snp", "beta", "se", "t_stat", "pvalue", "fdr",
                           "indep", "conditioning_snps")]
    }
    
    allbutone_tables[[length(allbutone_tables) + 1]] <- result
  }
  
  return(allbutone_tables)
}


# =========================================================================
# HELPER FUNCTION: Test SNPs with indices
# =========================================================================

#' Test SNPs with given column indices using univariate linear regression
#'
#' @param bigsnp bigSNP object with \code{$genotypes} (FBM) and \code{$map}
#' @param y Numeric phenotype vector (length = length(ind.row))
#' @param snp_indices Integer vector of column indices into the genotype FBM
#' @param snp_names Character vector of SNP names (same length as snp_indices)
#' @param design_base Data frame of covariates (row names = sample IDs)
#' @param ind.row Integer vector of row indices for samples in the genotype FBM
#' @param snp_conditioning Character vector of SNP names to condition on (default NULL)
#' @param ncores Number of cores for parallel computation (default 1)
#'
#' @return Data frame with columns: snp, beta, se, t_stat, pvalue, fdr
#' @export
test_snps_with_indices <- function(bigsnp, y, snp_indices, snp_names,
                                   design_base, ind.row,
                                   snp_conditioning = NULL,
                                   ncores = 1) {
  
  # Augment covariates with conditioning SNPs if provided
  if (!is.null(snp_conditioning) && length(snp_conditioning) > 0) {
    covar_df <- add_snps_to_covariates(
      bigsnp    = bigsnp,
      snp_names = snp_conditioning,
      covar_df  = design_base
    )
  } else {
    covar_df <- design_base
  }
  
  fit <- bigstatsr::big_univLinReg(
    X           = bigsnp$genotypes,
    y.train     = y,
    ind.train   = ind.row,
    ind.col     = snp_indices,
    covar.train = bigstatsr::covar_from_df(covar_df),
    ncores      = ncores
  )
  
  pvals <- stats::predict(fit, log10 = FALSE)
  BHfdr <- stats::p.adjust(pvals, "BH")
  
  data.frame(
    snp    = snp_names,
    beta   = fit$estim,
    se     = fit$std.err,
    t_stat = fit$score,
    pvalue = pvals,
    fdr    = BHfdr,
    stringsAsFactors = FALSE
  )
}