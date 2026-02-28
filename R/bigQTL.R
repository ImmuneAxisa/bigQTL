# =========================================================================
# MAIN FUNCTION: Conditional eQTL analysis
# =========================================================================

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
  
  library(bigstatsr)
  library(tidyverse)
  library(arrow)
  
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
  
  stepwise_dataset <- open_dataset(stepwise_dir)
  
  allbutone_dataset <- if (dir.exists(allbutone_dir) && length(list.files(allbutone_dir)) > 0) {
    open_dataset(allbutone_dir)
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
  ) %>%
    mutate(step = 0, conditioning_snps = NA_character_, .before = everything())
  
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
  write_parquet(bind_rows(stepwise_all), file.path(sw_dir, "part-0.parquet"))
  message(sprintf("    Wrote stepwise results to %s", gene_partition))
  
  # All-but-one
  if (length(allbutone_all) > 0) {
    abo_dir <- file.path(allbutone_dir, gene_partition)
    dir.create(abo_dir, recursive = TRUE, showWarnings = FALSE)
    write_parquet(bind_rows(allbutone_all), file.path(abo_dir, "part-0.parquet"))
    message(sprintf("    Wrote all-but-one results to %s", gene_partition))
  }
  
  return(invisible(NULL))
}


# =========================================================================
# STEPWISE CONDITIONING
# =========================================================================

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
  lead <- results_step0 %>% arrange(pvalue) %>% slice(1)
  
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
    ) %>%
      mutate(
        step = step,
        conditioning_snps = paste(conditioning_snps, collapse = ";"),
        .before = everything()
      )
    
    stepwise_tables[[length(stepwise_tables) + 1]] <- results_step
    
    new_lead <- results_step %>% arrange(pvalue) %>% slice(1)
    
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
      result <- stepwise_tables[[length(stepwise_tables)]] %>%
        mutate(
          indep = i,
          conditioning_snps = paste(snps_condition_on, collapse = ";")
        ) %>%
        select(snp, beta, se, t_stat, pvalue, fdr, indep, conditioning_snps)
      
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
      ) %>%
        mutate(
          indep = i,
          conditioning_snps = paste(snps_condition_on, collapse = ";"),
          .before = "snp"
        ) %>%
        select(snp, beta, se, t_stat, pvalue, fdr, indep, conditioning_snps)
    }
    
    allbutone_tables[[length(allbutone_tables) + 1]] <- result
  }
  
  return(allbutone_tables)
}


# =========================================================================
# HELPER FUNCTION: Test SNPs with indices
# =========================================================================

test_snps_with_indices <- function(bigsnp, y, snp_indices, snp_names,
                                   design_base, ind.row,
                                   snp_conditioning = NULL,
                                   ncores = 1) {
  
  library(bigstatsr)
  
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
  
  fit <- big_univLinReg(
    X           = bigsnp$genotypes,
    y.train     = y,
    ind.train   = ind.row,
    ind.col     = snp_indices,
    covar.train = covar_from_df(covar_df),
    ncores      = ncores
  )
  
  pvals <- predict(fit, log10 = FALSE)
  BHfdr <- p.adjust(pvals, "BH")
  
  tibble(
    snp    = snp_names,
    beta   = fit$estim,
    se     = fit$std.err,
    t_stat = fit$score,
    pvalue = pvals,
    fdr    = BHfdr
  )
}