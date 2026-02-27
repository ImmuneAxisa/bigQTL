# =========================================================================
# MAIN FUNCTION: Conditional eQTL analysis
# =========================================================================

run_conditional_eqtl <- function(
    bigfeatures,               # bigFeatures object with expression data (samples x features)
    bigsnp,                    # bigSNP object with genotype data (samples x snps)
    features_coord,            # Data frame with feature_name, chromosome, start, end columns
    design_base,               # Design matrix with covariates (samples x covariates)
    ind.row = NULL,            # Integer vector of row indices to use (NULL = all)
    cis_window = 1e6,          # Padding around feature coordinates
    do_conditioning = TRUE,
    pval_threshold = 1e-3,
    do_allbutone = TRUE,
    do_rint = TRUE,            # Apply RINT to phenotype per gene (default TRUE)
    ncores = 1,
    output_dir = "./eqtl_results") {
  
  library(bigstatsr)
  library(tidyverse)
  library(arrow)
  if (do_rint) library(RNOmni)
  
  # If no ind.row provided, use all rows
  if (is.null(ind.row)) {
    ind.row <- rows_along(bigsnp$genotypes)
  }
  
  # Validate features_coord
  required_cols <- c("feature_name", "chromosome", "start", "end")
  if (!all(required_cols %in% colnames(features_coord))) {
    stop(sprintf("features_coord must contain columns: %s",
                 paste(required_cols, collapse = ", ")))
  }
  
  # Create output directories
  stepwise_dir <- file.path(output_dir, "stepwise")
  allbutone_dir <- file.path(output_dir, "allbutone")
  
  if (!dir.exists(stepwise_dir)) {
    dir.create(stepwise_dir, recursive = TRUE)
  }
  if (!dir.exists(allbutone_dir)) {
    dir.create(allbutone_dir, recursive = TRUE)
  }
  
  genes <- features_coord$feature_name
  
  for (gene in genes) {
    
    message(sprintf("Processing %s...", gene))
    
    # Get gene index and coordinates
    gene_row <- features_coord[features_coord$feature_name == gene, ]
    gene_idx <- get_feature_indices(bigfeatures, gene)
    
    # Extract phenotype for kept individuals only
    y <- bigfeatures$features[ind.row, gene_idx]
    
    # RINT transform
    if (do_rint) {
      y <- RankNorm(y)
    }
    
    gene_chr <- gene_row$chromosome
    gene_start <- gene_row$start
    gene_end <- gene_row$end
    
    # Get cis SNPs on the fly
    cis_result <- get_cis_snps(bigsnp, gene_chr, gene_start, gene_end, cis_window)
    snp_indices <- cis_result$indices
    cis_snps_gene <- cis_result$names
    
    if (length(snp_indices) == 0) {
      message(sprintf("  No cis SNPs found for %s", gene))
      next
    }
    
    message(sprintf("  Found %d cis SNPs", length(snp_indices)))
    
    # Initialize lists for this gene's results
    stepwise_gene <- list()
    allbutone_gene <- list()
    
    # =====================================================================
    # STEP 0: MARGINAL TESTING
    # =====================================================================
    
    results_step0 <- test_snps_with_indices(
      bigsnp = bigsnp,
      y = y,
      snp_indices = snp_indices,
      snp_names = cis_snps_gene,
      design_base = design_base,
      ind.row = ind.row,
      ncores = ncores
    )
    
    results_step0 <- results_step0 %>%
      mutate(
        step = 0,
        conditioning_snps = NA_character_,
        .before = everything()
      )
    
    stepwise_gene[[1]] <- results_step0
    
    # =====================================================================
    # STEPWISE CONDITIONING
    # =====================================================================
    
    conditioning_snps <- character()
    
    if (do_conditioning) {
      
      lead_snp_step0 <- results_step0 %>%
        arrange(pvalue) %>%
        slice(1)
      
      lead_snp_name <- lead_snp_step0$snp
      lead_snp_pval <- lead_snp_step0$pvalue
      
      if (lead_snp_pval < pval_threshold) {
        
        message(sprintf("  Step 0: Lead SNP %s passes (p = %.2e)",
                        lead_snp_name, lead_snp_pval))
        
        conditioning_snps <- lead_snp_name
        step <- 1
        
        while (TRUE) {
          
          # Build covariate matrix using helper function
          covar_with_conditioning <- add_snps_to_covariates(
            bigsnp = bigsnp,
            snp_names = conditioning_snps,
            covar_df = design_base,
            ind.row = ind.row
          )
          
          results_step <- test_snps_with_indices(
            bigsnp = bigsnp,
            y = y,
            snp_indices = snp_indices,
            snp_names = cis_snps_gene,
            design_base = covar_with_conditioning,
            ind.row = ind.row,
            ncores = ncores
          )
          
          results_step <- results_step %>%
            mutate(
              step = step,
              conditioning_snps = paste(conditioning_snps, collapse = ";"),
              .before = everything()
            )
          
          stepwise_gene[[length(stepwise_gene) + 1]] <- results_step
          
          lead_snp_step <- results_step %>%
            arrange(pvalue) %>%
            slice(1)
          
          new_lead_snp <- lead_snp_step$snp
          new_lead_pval <- lead_snp_step$pvalue
          
          if (new_lead_pval < pval_threshold) {
            
            message(sprintf("    Step %d: New lead SNP %s passes (p = %.2e)",
                            step, new_lead_snp, new_lead_pval))
            
            conditioning_snps <- c(conditioning_snps, new_lead_snp)
            step <- step + 1
            
          } else {
            
            message(sprintf("    Step %d: Stopping (min p = %.2e)",
                            step, new_lead_pval))
            break
          }
        }
      }
    }
    
    # =====================================================================
    # ALL-BUT-ONE CONDITIONING
    # =====================================================================
    
    if (do_allbutone & length(conditioning_snps) > 1) {
      
      message(sprintf("  Running all-but-one conditioning for %d SNPs",
                      length(conditioning_snps)))
      
      for (i in seq_along(conditioning_snps)) {
        
        snp_focal <- conditioning_snps[i]
        snps_condition_on <- conditioning_snps[-i]
        
        if (i == length(conditioning_snps)) {
          
          # Copy from final stepwise step
          final_stepwise_step <- stepwise_gene[[length(stepwise_gene)]]
          
          result_allbutone <- final_stepwise_step %>%
            mutate(
              indep = i,
              conditioning_snps = paste(snps_condition_on, collapse = ";")
            ) %>%
            select(snp, beta, se, t_stat, pvalue, fdr, indep, conditioning_snps)
          
        } else {
          
          # Build covariate matrix using helper function
          covar_allbutone <- add_snps_to_covariates(
            bigsnp = bigsnp,
            snp_names = snps_condition_on,
            covar_df = design_base,
            ind.row = ind.row
          )
          
          result_allbutone <- test_snps_with_indices(
            bigsnp = bigsnp,
            y = y,
            snp_indices = snp_indices,
            snp_names = cis_snps_gene,
            design_base = covar_allbutone,
            ind.row = ind.row,
            ncores = ncores
          )
          
          result_allbutone <- result_allbutone %>%
            mutate(
              indep = i,
              conditioning_snps = paste(snps_condition_on, collapse = ";"),
              .before = "snp"
            ) %>%
            select(snp, beta, se, t_stat, pvalue, fdr, indep, conditioning_snps)
        }
        
        allbutone_gene[[length(allbutone_gene) + 1]] <- result_allbutone
      }
    }
    
    # =====================================================================
    # WRITE THIS GENE'S RESULTS AS PARTITION
    # =====================================================================
    
    # Write stepwise as hive partition
    gene_partition_name <- paste0("gene=", gene)
    stepwise_gene_dir <- file.path(stepwise_dir, gene_partition_name)
    dir.create(stepwise_gene_dir, recursive = TRUE)
    
    write_parquet(
      bind_rows(stepwise_gene),
      file.path(stepwise_gene_dir, "part-0.parquet")
    )
    
    message(sprintf("    Wrote stepwise results to %s", gene_partition_name))
    
    # Write all-but-one if it exists
    if (length(allbutone_gene) > 0) {
      
      allbutone_gene_dir <- file.path(allbutone_dir, gene_partition_name)
      dir.create(allbutone_gene_dir, recursive = TRUE)
      
      write_parquet(
        bind_rows(allbutone_gene),
        file.path(allbutone_gene_dir, "part-0.parquet")
      )
      
      message(sprintf("    Wrote all-but-one results to %s", gene_partition_name))
    }
  }
  
  # =====================================================================
  # OPEN AND RETURN PARQUET DATASETS
  # =====================================================================
  
  message("\nOpening Parquet datasets...")
  
  # Stepwise dataset
  stepwise_dataset <- open_dataset(stepwise_dir)
  
  # All-but-one dataset
  allbutone_dataset <- if (dir.exists(allbutone_dir) && length(list.files(allbutone_dir)) > 0) {
    open_dataset(allbutone_dir)
  } else {
    NULL
  }
  
  message(sprintf("Results saved to: %s", output_dir))
  
  return(list(
    stepwise = stepwise_dataset,
    allbutone = allbutone_dataset,
    stepwise_dir = stepwise_dir,
    allbutone_dir = allbutone_dir
  ))
}


# =========================================================================
# HELPER FUNCTION: Test SNPs with indices
# =========================================================================

test_snps_with_indices <- function(bigsnp, y, snp_indices, snp_names,
                                   design_base, ind.row = NULL,
                                   ncores = 1) {
  
  library(bigstatsr)
  
  if (is.null(ind.row)) {
    ind.row <- rows_along(bigsnp$genotypes)
  }
  
  fit <- big_univLinReg(
    X           = bigsnp$genotypes,
    y.train     = y,
    ind.train   = ind.row,
    ind.col     = snp_indices,
    covar.train = covar_from_df(design_base),
    ncores      = ncores
  )
  
  pvals <- predict(fit, log10 = FALSE)
  BHfdr <- p.adjust(pvals, "BH")
  
  results <- tibble(
    snp = snp_names,
    beta = fit$estim,
    se = fit$std.err,
    t_stat = fit$score,
    pvalue = pvals,
    fdr = BHfdr
  )
  
  return(results)
}

