# GitHub Copilot Instructions for bigQTL

## Package Overview

**bigQTL** is an R package for scalable conditional cis-QTL analysis using
file-backed matrices from the `bigstatsr`/`bigsnpr` ecosystem. Key capabilities:

- **Conditional eQTL analysis** – stepwise conditioning and all-but-one
  conditioning via `run_conditional_eqtl()`
- **eigenMT multiple-testing correction** – per-gene effective number of
  independent tests using Ledoit-Wolf shrinkage (`eigenMT_gene()`,
  `eigenMT_batch()`)
- **Data-preparation utilities** – genotype PCA (`compute_geno_pcs()`),
  feature PCA (`compute_feature_pcs()`), RINT transformation (`rint()`)
- **bigFeatures S3 class** – wraps expression/feature matrices in an FBM
  for disk-backed storage

Results are written as partitioned Apache Parquet datasets (via `arrow`).

## Repository Layout

```
R/
  bigQTL.R          # run_conditional_eqtl(), process_gene(), run_stepwise(),
                    #   run_allbutone(), test_snps_with_indices()
  bigFeatures.R     # bigFeatures() constructor
  helpers.R         # add_snps_to_covariates(), get_cis_snps(),
                    #   get_feature_indices()
  eigenMT.R         # lw_shrink_cor(), count_eigenvalues(), eigenMT_correct(),
                    #   eigenMT_gene(), eigenMT_batch()
  prep_helpers.R    # compute_geno_pcs(), compute_feature_pcs(), rint()
  bigQTL-package.R  # package-level documentation
vignettes/
  articles/
    eigenMT_validation.Rmd  # interactive validation of R vs Python eigenMT
DESCRIPTION
NAMESPACE
```

## Environment Setup

### System requirements

- **R ≥ 4.1**
- Suggested system libraries: `libhdf5-dev`, `libcurl4-openssl-dev`

### Install R dependencies

```r
# Install package dependencies from DESCRIPTION
install.packages(c("bigstatsr", "bigsnpr", "arrow"))

# Suggested (needed for eigenMT validation vignette)
install.packages(c("reticulate", "knitr", "rmarkdown"))
```

### Python environment (eigenMT validation only)

The `vignettes/articles/eigenMT_validation.Rmd` vignette compares the R
eigenMT implementation against the original Python version. Set up the conda
environment once:

```r
library(reticulate)
reticulate::conda_create(
  "eigenMT",
  packages = c("numpy", "scipy", "scikit-learn", "pandas"),
  python_version = "3.8"
)
```

## Building & Checking the Package

```bash
# From the repository root (shell)
R CMD build .
R CMD check bigQTL_*.tar.gz
```

```r
# From an R session
devtools::document()   # regenerate NAMESPACE and Rd files
devtools::build()
devtools::check()
```

## Running Tests

There is currently no `tests/testthat/` suite. Functional validation lives in
`vignettes/articles/eigenMT_validation.Rmd` (all chunks use `eval = FALSE` and
are intended for interactive execution).

To run the inline unit-test section of that vignette interactively:

```r
library(bigQTL)
library(bigstatsr)
# ... follow the vignette setup steps, then run the "Unit Tests" chunk
```

The CI workflow (`.github/workflows/R-CMD-check.yaml`) runs `R CMD check`
on every push and pull request, which exercises `devtools::check()` including
example code in Rd files and vignette build (skipped if `eval = FALSE`).

## Key Conventions

- **File-backed matrices (FBM)**: All genotype and feature data live on disk
  as `bigstatsr::FBM` objects; never materialize large matrices into R memory.
- **Sample ID matching**: `design_base` row names must be meaningful sample IDs
  that appear in both `bigsnp$fam$sample.ID` and
  `bigfeatures$rowData$sample_name`.
- **Output format**: Results are written as partitioned Parquet under
  `output_dir/stepwise/gene=<name>/part-0.parquet` and
  `output_dir/allbutone/gene=<name>/part-0.parquet`.
- **Parallelism**: Within-gene SNP regression uses `ncores`; across-gene
  parallelism uses `ncores_genes` (both default to 1; increase on
  multi-core machines).
- **RINT**: By default (`do_rint = TRUE`) the phenotype is rank inverse-normal
  transformed per gene before regression.

## Adding New Functions

1. Place source in the relevant `R/*.R` file (or a new file if the scope
   warrants it).
2. Add roxygen2 documentation with `@param`, `@return`, and `@export`
   (if public) or `@keywords internal` (if private).
3. Run `devtools::document()` to update `NAMESPACE` and Rd files.
4. Add an example to the Rd block or extend the vignette.
