#!/usr/bin/env Rscript
# =========================================================================
# Generate a realistic synthetic QTL dataset with known ground truth.
#
# 5 phenotypes on the same chromosome, each with its own non-overlapping
# cis-window, carrying 3, 2, 1, 0, 0 independent causal signals
# respectively. Each causal signal has a couple of "LD friend" SNPs
# (variable r^2 in [0.2, 0.9]) that tag it imperfectly, the way a real
# genotyped SNP tags an untyped causal variant -- these are what stepwise
# conditioning must learn to discard once the true causal SNP is in the
# model.
#
# Output:
#   - R/sysdata.rda: internal `dummy_qtl_data` list (raw matrices + map +
#     ground truth, no FBMs -- FBM backing files don't survive
#     serialization, see tests/testthat/helper-dummy-data.R for how the
#     bigsnp/bigpheno objects are reconstructed from this at test time).
#   - dev/dummy_data_association_plots.png: one panel per phenotype,
#     one line per stepwise step (i.e. per recovered signal), for visual
#     inspection that the simulated signals are recoverable.
#
# Re-run this script to regenerate both outputs; it is fully seeded.
# =========================================================================

devtools::load_all(".", quiet = TRUE)
library(ggplot2)

SEED <- 20260929
set.seed(SEED)

# =========================================================================
# Simulation parameters
# =========================================================================

CHROM                <- "1"
N_SAMPLES            <- 400
N_SIGNALS            <- c(3, 2, 1, 0, 0)          # per phenotype
N_FRIENDS_PER_SIGNAL <- 2
N_FILLER_SNPS        <- 60                        # per window
WINDOW_SPACING       <- 3e6                       # gap between gene bodies
GENE_WIDTH           <- 1e4                       # width of gene body
FLANK                <- 4e5                       # SNPs within +/- FLANK of
                                                   # gene body (< default
                                                   # 1e6 cis_window, and
                                                   # << WINDOW_SPACING/2 so
                                                   # windows never overlap)
VAR_EXPLAINED_RANGE  <- c(0.08, 0.18)              # per-causal-SNP variance
FRIEND_R2_RANGE      <- c(0.2, 0.9)
CIS_WINDOW           <- 1e6                        # matches package default

# Background "genome-wide backbone" SNPs and non-tested phenotypes: these
# carry no simulated genetic effect and are never cis-SNPs for any of the
# 5 test phenotypes (placed far away / excluded from pheno_coord). Their
# only purpose is to keep bigQTL()'s compute_geno_pcs()/compute_pheno_pcs()
# realistic -- without them the genotype/phenotype PCs would be computed
# from nothing but the 5 causal test loci/traits themselves and would
# partially soak up the very signal we're trying to detect (real datasets
# compute PCs from thousands of genes / hundreds of thousands of SNPs, so
# a handful of causal loci/traits are a negligible fraction of the basis).
N_BACKGROUND_SNPS    <- 3000
N_BACKGROUND_PHENOS  <- 200
BACKGROUND_SNP_START <- 5e7   # far beyond every phenotype's cis window

n_pheno <- length(N_SIGNALS)
sample_ids <- sprintf("sample_%03d", seq_len(N_SAMPLES))

# =========================================================================
# Genotype simulation helpers
# =========================================================================

# Simulate a biallelic SNP under Hardy-Weinberg equilibrium at a given MAF.
# Returns the two haplotype-allele vectors plus the dosage (0/1/2).
sim_hap_snp <- function(n, maf) {
  h1 <- stats::rbinom(n, 1, maf)
  h2 <- stats::rbinom(n, 1, maf)
  list(h1 = h1, h2 = h2, dosage = h1 + h2)
}

# Simulate an "LD friend" of a causal SNP via a haplotype-copying model:
# each haplotype is copied from the causal SNP's haplotype with probability
# c = sqrt(target_r2), and redrawn independently (at the same MAF)
# otherwise. Because the two haplotypes are copied independently and the
# marginal allele frequency is preserved, cor(friend_dosage, causal_dosage)
# equals c in expectation, i.e. the realized r^2 targets target_r2.
sim_ld_friend <- function(causal_h1, causal_h2, maf, target_r2) {
  n <- length(causal_h1)
  c_copy <- sqrt(target_r2)
  copy1 <- stats::rbinom(n, 1, c_copy) == 1
  copy2 <- stats::rbinom(n, 1, c_copy) == 1
  f1 <- ifelse(copy1, causal_h1, stats::rbinom(n, 1, maf))
  f2 <- ifelse(copy2, causal_h2, stats::rbinom(n, 1, maf))
  f1 + f2
}

# =========================================================================
# Simulate one phenotype's cis-window: filler SNPs + causal SNPs + their
# LD friends. Returns the genotype columns, map rows, and ground-truth
# records for this window.
# =========================================================================

simulate_window <- function(pheno_name, gene_start, gene_end, n_signals,
                            n_samples, sample_ids) {

  window_lo <- gene_start - FLANK
  window_hi <- gene_end + FLANK

  n_causal  <- n_signals
  n_friends <- n_causal * N_FRIENDS_PER_SIGNAL
  n_total   <- N_FILLER_SNPS + n_causal + n_friends

  # Reserve unique integer positions for every SNP in this window up front
  # so filler/causal/friend positions never collide.
  all_positions <- sample(window_lo:window_hi, n_total, replace = FALSE)
  pos_causal  <- if (n_causal > 0) sort(all_positions[seq_len(n_causal)]) else integer(0)
  pos_friends <- if (n_friends > 0) all_positions[n_causal + seq_len(n_friends)] else integer(0)
  pos_filler  <- all_positions[n_causal + n_friends + seq_len(N_FILLER_SNPS)]

  geno_cols  <- list()
  map_rows   <- list()
  causal_gt  <- list()
  friend_gt  <- list()

  # ---- filler SNPs (no effect on phenotype) ----
  for (j in seq_along(pos_filler)) {
    maf <- stats::runif(1, 0.05, 0.5)
    snp <- sim_hap_snp(n_samples, maf)
    id  <- sprintf("chr%s_%d_filler", CHROM, pos_filler[j])
    geno_cols[[id]] <- snp$dosage
    map_rows[[id]]  <- data.frame(chromosome = CHROM, marker.ID = id,
                                  physical.pos = pos_filler[j],
                                  stringsAsFactors = FALSE)
  }

  # ---- causal SNPs + their LD friends ----
  friend_counter <- 0L
  for (k in seq_len(n_causal)) {
    maf <- stats::runif(1, 0.15, 0.45)
    causal <- sim_hap_snp(n_samples, maf)
    id <- sprintf("chr%s_%d_causal_%s_%d", CHROM, pos_causal[k], pheno_name, k)
    geno_cols[[id]] <- causal$dosage
    map_rows[[id]]  <- data.frame(chromosome = CHROM, marker.ID = id,
                                  physical.pos = pos_causal[k],
                                  stringsAsFactors = FALSE)

    var_explained <- stats::runif(1, VAR_EXPLAINED_RANGE[1], VAR_EXPLAINED_RANGE[2])
    beta <- sqrt(var_explained / stats::var(causal$dosage)) * sample(c(-1, 1), 1)

    causal_gt[[id]] <- data.frame(
      pheno_name = pheno_name, snp = id, maf = maf, beta = beta,
      var_explained = var_explained, stringsAsFactors = FALSE
    )

    for (f in seq_len(N_FRIENDS_PER_SIGNAL)) {
      friend_counter <- friend_counter + 1L
      pos_f <- pos_friends[friend_counter]
      target_r2 <- stats::runif(1, FRIEND_R2_RANGE[1], FRIEND_R2_RANGE[2])
      friend_dosage <- sim_ld_friend(causal$h1, causal$h2, maf, target_r2)
      fid <- sprintf("chr%s_%d_friend_%s_%d_%d", CHROM, pos_f, pheno_name, k, f)
      geno_cols[[fid]] <- friend_dosage
      map_rows[[fid]]  <- data.frame(chromosome = CHROM, marker.ID = fid,
                                     physical.pos = pos_f,
                                     stringsAsFactors = FALSE)
      achieved_r2 <- stats::cor(friend_dosage, causal$dosage)^2
      friend_gt[[fid]] <- data.frame(
        pheno_name = pheno_name, snp = fid, causal_snp = id,
        target_r2 = target_r2, achieved_r2 = achieved_r2,
        stringsAsFactors = FALSE
      )
    }
  }

  geno_mat <- do.call(cbind, geno_cols)
  rownames(geno_mat) <- sample_ids
  map_df <- do.call(rbind, map_rows)
  map_df <- map_df[order(map_df$physical.pos), ]
  geno_mat <- geno_mat[, map_df$marker.ID, drop = FALSE]

  list(
    geno_mat = geno_mat,
    map      = map_df,
    causal   = if (length(causal_gt) > 0) do.call(rbind, causal_gt) else NULL,
    friends  = if (length(friend_gt) > 0) do.call(rbind, friend_gt) else NULL
  )
}

# =========================================================================
# Build all 5 windows + phenotype coordinates
# =========================================================================

windows       <- vector("list", n_pheno)
pheno_names   <- sprintf("pheno%d", seq_len(n_pheno))
pheno_coord_rows <- vector("list", n_pheno)

for (i in seq_len(n_pheno)) {
  gene_start <- 1e6 + (i - 1) * WINDOW_SPACING
  gene_end   <- gene_start + GENE_WIDTH

  windows[[i]] <- simulate_window(
    pheno_name = pheno_names[i], gene_start = gene_start, gene_end = gene_end,
    n_signals = N_SIGNALS[i], n_samples = N_SAMPLES, sample_ids = sample_ids
  )

  pheno_coord_rows[[i]] <- data.frame(
    pheno_name = pheno_names[i], chromosome = CHROM,
    start = gene_start, end = gene_end, stringsAsFactors = FALSE
  )
}

geno_mat <- do.call(cbind, lapply(windows, `[[`, "geno_mat"))
map_df   <- do.call(rbind, lapply(windows, `[[`, "map"))
pheno_coord <- do.call(rbind, pheno_coord_rows)

# =========================================================================
# Background genotype backbone (no effect on any phenotype; see param
# comments above). Independent SNPs, no LD structure needed.
# =========================================================================

bg_positions <- sort(sample(BACKGROUND_SNP_START + seq_len(N_BACKGROUND_SNPS * 100),
                            N_BACKGROUND_SNPS, replace = FALSE))
bg_geno_cols <- lapply(bg_positions, function(pos) {
  sim_hap_snp(N_SAMPLES, stats::runif(1, 0.05, 0.5))$dosage
})
bg_ids <- sprintf("chr%s_%.0f_background", CHROM, bg_positions)
bg_geno_mat <- do.call(cbind, bg_geno_cols)
dimnames(bg_geno_mat) <- list(sample_ids, bg_ids)
bg_map_df <- data.frame(chromosome = CHROM, marker.ID = bg_ids,
                        physical.pos = bg_positions, stringsAsFactors = FALSE)

geno_mat <- cbind(geno_mat, bg_geno_mat)
map_df   <- rbind(map_df, bg_map_df)
map_df   <- map_df[order(map_df$physical.pos), ]
geno_mat <- geno_mat[, map_df$marker.ID, drop = FALSE]

# =========================================================================
# Simulate phenotype values from each window's causal SNPs + noise
# =========================================================================

pheno_mat <- matrix(NA_real_, nrow = N_SAMPLES, ncol = n_pheno,
                    dimnames = list(sample_ids, pheno_names))

for (i in seq_len(n_pheno)) {
  causal <- windows[[i]]$causal
  if (is.null(causal)) {
    y_signal <- 0
  } else {
    causal_geno <- windows[[i]]$geno_mat[, causal$snp, drop = FALSE]
    y_signal <- as.vector(causal_geno %*% causal$beta)
  }
  total_var_explained <- if (is.null(causal)) 0 else sum(causal$var_explained)
  noise_sd <- sqrt(max(0.05, 1 - total_var_explained))
  pheno_mat[, i] <- y_signal + stats::rnorm(N_SAMPLES, sd = noise_sd)
}

# ---- background phenotypes (pure noise, not in pheno_coord -- never
# cis-QTL-tested, only there to keep compute_pheno_pcs() realistic) ----

bg_pheno_names <- sprintf("bg_pheno_%03d", seq_len(N_BACKGROUND_PHENOS))
bg_pheno_mat <- matrix(stats::rnorm(N_SAMPLES * N_BACKGROUND_PHENOS),
                       nrow = N_SAMPLES, dimnames = list(sample_ids, bg_pheno_names))
pheno_mat <- cbind(pheno_mat, bg_pheno_mat)

# =========================================================================
# Design matrix: a couple of covariates uncorrelated with genotype/phenotype
# =========================================================================

design_base <- data.frame(
  age = round(stats::rnorm(N_SAMPLES, mean = 50, sd = 10), 1),
  sex = factor(sample(c("F", "M"), N_SAMPLES, replace = TRUE)),
  row.names = sample_ids
)

# =========================================================================
# Ground truth
# =========================================================================

ground_truth <- list(
  n_signals            = stats::setNames(N_SIGNALS, pheno_names),
  causal               = do.call(rbind, lapply(windows, `[[`, "causal")),
  friends              = do.call(rbind, lapply(windows, `[[`, "friends")),
  n_background_snps    = N_BACKGROUND_SNPS,
  n_background_phenos  = N_BACKGROUND_PHENOS
)
rownames(ground_truth$causal)  <- NULL
rownames(ground_truth$friends) <- NULL

cat(sprintf("Simulated %d SNPs across %d phenotypes (%d samples)\n",
           nrow(map_df), n_pheno, N_SAMPLES))
cat("Causal SNPs:\n")
print(ground_truth$causal[, c("pheno_name", "snp", "maf", "beta", "var_explained")])
cat("\nLD friends (target vs. achieved r^2):\n")
print(ground_truth$friends[, c("pheno_name", "snp", "causal_snp", "target_r2", "achieved_r2")])

# =========================================================================
# Sanity check: recover signals with the package's own stepwise conditioning
# and generate one association-profile line per phenotype/step for
# inspection. This also doubles as a smoke test that the ground truth is
# actually recoverable before we freeze it as a test fixture.
# =========================================================================

geno_fbm <- bigstatsr::as_FBM(geno_mat)
bigsnp_check <- list(
  genotypes = geno_fbm,
  fam       = data.frame(sample.ID = sample_ids, stringsAsFactors = FALSE),
  map       = map_df
)
ind_row <- seq_len(N_SAMPLES)

plot_data <- list()
recovered_summary <- list()

for (i in seq_len(n_pheno)) {
  pheno_name <- pheno_names[i]
  y <- rint(pheno_mat[, i])

  cis <- get_cis_snps(bigsnp_check, pheno_chr = CHROM,
                      pheno_start = pheno_coord$start[i],
                      pheno_end = pheno_coord$end[i], cis_window = CIS_WINDOW)

  results_step0 <- test_snps_with_indices(
    bigsnp = bigsnp_check, y = y, snp_indices = cis$indices,
    snp_names = cis$names, design_base = design_base, ind.row = ind_row
  )
  results_step0$step <- 0L
  results_step0$conditioning_snps <- NA_character_

  stepwise <- run_stepwise(
    results_step0 = results_step0, bigsnp = bigsnp_check, y = y,
    snp_indices = cis$indices, cis_snps_pheno = cis$names,
    design_base = design_base, ind.row.snp = ind_row,
    pval_threshold = 1e-6, max_steps = 5, ncores = 1, pheno = pheno_name
  )

  recovered_summary[[pheno_name]] <- length(stepwise$conditioning_snps)

  for (s in seq_along(stepwise$stepwise_tables)) {
    tbl <- stepwise$stepwise_tables[[s]]
    tbl <- merge(tbl, map_df[, c("marker.ID", "physical.pos")],
                by.x = "snp", by.y = "marker.ID")
    tbl$pheno_name <- pheno_name
    tbl$step_label <- sprintf("step %d", tbl$step[1])
    plot_data[[length(plot_data) + 1]] <- tbl
  }
}

cat("\nRecovered independent signals (stepwise) vs. ground truth:\n")
print(data.frame(
  pheno_name = pheno_names,
  n_signals_truth = N_SIGNALS,
  n_signals_recovered = unlist(recovered_summary[pheno_names])
))

plot_df <- do.call(rbind, plot_data)
plot_df$pheno_name <- factor(plot_df$pheno_name, levels = pheno_names)

causal_marks <- ground_truth$causal
causal_marks <- merge(causal_marks, map_df[, c("marker.ID", "physical.pos")],
                      by.x = "snp", by.y = "marker.ID")
causal_marks$pheno_name <- factor(causal_marks$pheno_name, levels = pheno_names)

p <- ggplot(plot_df, aes(x = physical.pos, y = -log10(pvalue), color = step_label)) +
  geom_line() +
  geom_vline(data = causal_marks, aes(xintercept = physical.pos),
            linetype = "dashed", color = "grey40") +
  facet_wrap(~ pheno_name, scales = "free_x", ncol = 1) +
  labs(x = "Position (bp)", y = expression(-log[10](p)), color = "Conditioning step",
      title = "Dummy QTL dataset: association profile per stepwise-conditioning step",
      subtitle = "Dashed lines mark true causal SNP positions") +
  theme_bw()

out_png <- file.path("dev", "dummy_data_association_plots.png")
ggsave(out_png, p, width = 9, height = 14, dpi = 150)
cat(sprintf("\nWrote association plots to %s\n", out_png))

# =========================================================================
# Full-pipeline sanity check: run bigQTL() itself (PC computation +
# augmentation + conditional QTL analysis together) and confirm the
# background scaffold above is large enough that adding genotype/phenotype
# PCs as covariates doesn't eat into the simulated signal.
# =========================================================================

geno_fbm_code <- bigstatsr::FBM.code256(
  nrow = nrow(geno_mat), ncol = ncol(geno_mat), code = bigsnpr::CODE_012
)
geno_fbm_code[, ] <- geno_mat
bigsnp_full <- list(
  genotypes = geno_fbm_code,
  fam       = data.frame(sample.ID = sample_ids, stringsAsFactors = FALSE),
  map       = map_df
)
bigpheno_full <- bigPheno(pheno_mat)

bigqtl_out_dir <- tempfile("dummy_qtl_bigQTL_check_")
# compute_pheno_pcs()/compute_geno_pcs() use randomized SVD algorithms
# (big_randomSVD/snp_autoSVD); reseed immediately before the call so this
# check -- and any test calling bigQTL() on this fixture -- is
# reproducible. Tests must set.seed() the same way before calling bigQTL().
set.seed(1)
bigqtl_res <- bigQTL(
  bigpheno = bigpheno_full, bigsnp = bigsnp_full, pheno_coord = pheno_coord,
  design_base = design_base, n_pheno_pcs = 3, n_geno_pcs = 3,
  min_snps = 10, output_dir = bigqtl_out_dir, verbose = FALSE
)
bigqtl_df <- as.data.frame(bigqtl_res$stepwise)
# Number of recovered signals == length(conditioning_snps) inside
# run_stepwise(), which is NOT simply max(step): step 0 (marginal) is
# always written to the output even when its own lead SNP fails
# pval_threshold (i.e. zero real signals), so max(step) == 0 is ambiguous
# between "0 signals" and "exactly 1 signal, no further conditioning step
# passed". Any step >= 1 present unambiguously implies step 0 passed.
count_recovered <- function(df_pheno, pval_threshold = 1e-3) {
  max_step <- max(df_pheno$step)
  if (max_step >= 1) return(as.integer(max_step) + 1L)
  step0 <- df_pheno[df_pheno$step == 0, ]
  if (min(step0$pvalue) < pval_threshold) 1L else 0L
}
n_recovered_full <- vapply(pheno_names, function(p) {
  count_recovered(bigqtl_df[bigqtl_df$pheno == p, ])
}, integer(1))

cat("\nRecovered independent signals via full bigQTL() pipeline (incl. PCs):\n")
print(data.frame(pheno_name = pheno_names, n_signals_truth = N_SIGNALS,
                 n_signals_recovered = n_recovered_full))
unlink(bigqtl_out_dir, recursive = TRUE)

if (!identical(unname(n_recovered_full), as.integer(N_SIGNALS))) {
  warning("Full bigQTL() pipeline did not recover the expected number of ",
         "signals -- consider a larger background scaffold or stronger effects.")
}

# =========================================================================
# Save as an internal package dataset (no FBMs -- those don't survive
# serialization; see tests/testthat/helper-dummy-data.R for reconstruction)
# =========================================================================

dummy_qtl_data <- list(
  geno         = geno_mat,
  map          = map_df,
  fam          = data.frame(sample.ID = sample_ids, stringsAsFactors = FALSE),
  pheno        = pheno_mat,
  pheno_coord  = pheno_coord,
  design_base  = design_base,
  ground_truth = ground_truth,
  seed         = SEED
)

usethis::use_data(dummy_qtl_data, internal = TRUE, overwrite = TRUE)
cat("Saved dummy_qtl_data to R/sysdata.rda\n")
