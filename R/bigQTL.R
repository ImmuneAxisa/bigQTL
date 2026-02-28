# =========================================================================
# MAIN FUNCTION: Conditional eQTL analysis
# =========================================================================

#' Run conditional eQTL analysis
#'
#' Performs conditional cis-QTL analysis using stepwise conditioning and
#' all-but-one conditioning on file-backed genotype and feature matrices.
#'
#' @param bigfeatures bigFeatures object with expression data (samples x features)
#' @param bigsnp bigSNP object with genotype data (samples x snps)
#' @param features_coord Data frame with feature_name, chromosome, start, end columns
#' @param design_base Data frame with covariates; row names = sample IDs
#' @param cis_window Padding around feature coordinates (default 1e6)
#' @param do_conditioning Logical; perform stepwise conditioning (default TRUE)
#' @param pval_threshold P-value threshold for stepwise conditioning (default 1e-3)
#' @param do_allbutone Logical; perform all-but-one conditioning (default TRUE)
#' @param do_rint Logical; apply RINT to phenotype per gene (default TRUE)
#' @param ncores Cores for per-SNP regression within-gene (default 1)
#' @param ncores_genes Cores for across-gene parallelisation (default 1)
#' @param output_dir Directory for Parquet output (default "./eqtl_results")
#'
#' @return List with elements: stepwise (Arrow dataset), allbutone (Arrow dataset or NULL),
#'   stepwise_dir, allbutone_dir
#' @export
run_conditional_eqtl <- function(
    bigfeatures,               # bigFeatures object with expression data (samples x features)
    bigsnp,                    # bigSNP object with genotype data (samples x snps)
    features_coord,            # Data frame with feature_name, chromosome, start, end columns
    design_base,               # Data frame with covariates; row names = sample IDs
    cis_window = 1e6,          # Padding around feature coordinates
    do_conditioning = TRUE,
    pval_threshold = 1e-3,
    do_allbutone = TRUE,
    do_rint = TRUE,            # Apply RINT to phenotype per gene (default TRUE)
    ncores = 1,                # Cores for per-SNP regression (within-gene)
    ncores_genes = 1,          # Cores for across-gene parallelisation
    output_dir = "./eqtl_results") {
  
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
  
  # Map to feature row indices
  ind.row.feat <- match(sample_ids, bigfeatures$rowData$sample_name)
  if (any(is.na(ind.row.feat))) {
    missing <- sample_ids[is.na(ind.row.feat)]
    stop(sprintf("%d sample ID(s) from design_base not found in bigfeatures$rowData$sample_name: %s",
                 length(missing), paste(head(missing, 3), collapse = ", ")))
  }
  
  message(sprintf("  Matched %d samples from design_base to bigsnp and bigfeatures",
                  length(sample_ids)))
  
  # Validate features_coord
  required_cols <- c("feature_name", "chromosome", "start", "end")
  if (!all(required_cols %in% colnames(features_coord))) {
    stop(sprintf("features_coord must contain columns: %s",
                 paste(required_cols, collapse = ", ")))
  }
  
  # Create output directories
  stepwise_dir <- file.path(output_dir, "stepwise")
  allbutone_dir <- file.path(output_dir, "allbutone")
  
  if (!dir.exists(stepwise_dir)) dir.create(stepwise_dir, recursive = TRUE)
  if (!dir.exists(allbutone_dir)) dir.create(allbutone_dir, recursive = TRUE)
  
  # =====================================================================
  # Process genes (parallelised with mclapply or sequential with lapply)
  # =====================================================================
  
  genes <- features_coord$feature_name
  
  process_fn <- function(gene) {
    process_gene(
      gene           = gene,
      bigfeatures    = bigfeatures,
      bigsnp         = bigsnp,
      features_coord = features_coord,
      design_base    = design_base,
      ind.row.snp    = ind.row.snp,
      ind.row.feat   = ind.row.feat,
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
  
  if (ncores_genes > 1) {
    parallel::mclapply(genes, process_fn, mc.cores = ncores_genes)
  } else {
    lapply(genes, process_fn)
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
# PER-GENE FUNCTION: Process a single gene
# =========================================================================

#' Process a single gene for conditional eQTL analysis
#'
#' @param gene Gene name to process
#' @param bigfeatures bigFeatures object
#' @param bigsnp bigSNP object
#' @param features_coord Data frame with gene coordinates
#' @param design_base Data frame with covariates
#' @param ind.row.snp Row indices into genotype FBM for samples
#' @param ind.row.feat Row indices into feature FBM for samples
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
process_gene <- function(gene, bigfeatures, bigsnp, features_coord,
                         design_base, ind.row.snp, ind.row.feat,
                         cis_window, do_conditioning, pval_threshold,
                         do_allbutone, do_rint, ncores,
                         stepwise_dir, allbutone_dir) {
  
  message(sprintf("Processing %s...", gene))
  
  # Get gene coordinates and feature index
  gene_row <- features_coord[features_coord$feature_name == gene, ]
  gene_idx <- get_feature_indices(bigfeatures, gene)
  
  # Extract phenotype for kept individuals
  y <- bigfeatures$features[ind.row.feat, gene_idx]
  
  # RINT transform
  if (do_rint) {
    y <- rint(y)
  }
  
  # Get cis SNPs
  cis_result <- get_cis_snps(bigsnp, gene_row$chromosome,
                             gene_row$start, gene_row$end, cis_window)
  snp_indices <- cis_result$indices
  cis_snps_gene <- cis_result$names
  
  if (length(snp_indices) == 0) {
    message(sprintf("  No cis SNPs found for %s", gene))
    return(invisible(NULL))
  }
  
  message(sprintf("  Found %d cis SNPs", length(snp_indices)))
  
  # ==== Step 0: Marginal testing ====
  
  results_step0 <- test_snps_with_indices(
    bigsnp      = bigsnp,
    y           = y,
    snp_indices = snp_indices,
    snp_names   = cis_snps_gene,
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
    cis_snps_gene  = cis_snps_gene,
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
    cis_snps_gene     = cis_snps_gene,
    design_base       = design_base,
    ind.row.snp       = ind.row.snp,
    do_allbutone      = do_allbutone,
    ncores            = ncores
  )
  
  # ==== Write results ====
  
  gene_partition <- paste0("gene=", gene)
  
  # Stepwise
  sw_dir <- file.path(stepwise_dir, gene_partition)
  dir.create(sw_dir, recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(do.call(rbind, stepwise_all), file.path(sw_dir, "part-0.parquet"))
  message(sprintf("    Wrote stepwise results to %s", gene_partition))
  
  # All-but-one
  if (length(allbutone_all) > 0) {
    abo_dir <- file.path(allbutone_dir, gene_partition)
    dir.create(abo_dir, recursive = TRUE, showWarnings = FALSE)
    arrow::write_parquet(do.call(rbind, allbutone_all), file.path(abo_dir, "part-0.parquet"))
    message(sprintf("    Wrote all-but-one results to %s", gene_partition))
  }
  
  return(invisible(NULL))
}


# =========================================================================
# STEPWISE CONDITIONING
# =========================================================================

#' Run stepwise conditioning for a gene
#'
#' @param results_step0 Data frame of marginal association results (step 0)
#' @param bigsnp bigSNP object
#' @param y Numeric phenotype vector
#' @param snp_indices Integer vector of cis-SNP column indices
#' @param cis_snps_gene Character vector of cis-SNP names
#' @param design_base Data frame of covariates
#' @param ind.row.snp Integer vector of row indices for samples in genotype FBM
#' @param do_conditioning Logical; perform stepwise conditioning
#' @param pval_threshold P-value threshold for stepwise conditioning
#' @param ncores Number of cores
#'
#' @return List with elements: stepwise_tables, conditioning_snps
#' @keywords internal
run_stepwise <- function(results_step0, bigsnp, y, snp_indices,
                         cis_snps_gene, design_base, ind.row.snp,
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
      snp_names        = cis_snps_gene,
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

#' Run all-but-one conditioning for a gene
#'
#' @param conditioning_snps Character vector of independent SNP names
#' @param stepwise_tables List of data frames from stepwise conditioning
#' @param bigsnp bigSNP object
#' @param y Numeric phenotype vector
#' @param snp_indices Integer vector of cis-SNP column indices
#' @param cis_snps_gene Character vector of cis-SNP names
#' @param design_base Data frame of covariates
#' @param ind.row.snp Integer vector of row indices for samples in genotype FBM
#' @param do_allbutone Logical; perform all-but-one conditioning
#' @param ncores Number of cores
#'
#' @return List of data frames with all-but-one results
#' @keywords internal
run_allbutone <- function(conditioning_snps, stepwise_tables, bigsnp, y,
                          snp_indices, cis_snps_gene, design_base,
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
        snp_names        = cis_snps_gene,
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