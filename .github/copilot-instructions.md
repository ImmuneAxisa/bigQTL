# GitHub Copilot Instructions for bigQTL

## Package Overview

**bigQTL** is an R package for scalable conditional cis-QTL analysis built to operate with file-backed matrices (FBM) from the bigstatsr / bigsnpr ecosystem. Key capabilities:

- Conditional eQTL analysis — stepwise conditioning and all-but-one conditioning via `run_conditional_eqtl()`
- eigenMT multiple-testing correction — per-gene effective number of independent tests with Ledoit–Wolf shrinkage (`lw_shrink_geno()`, `eigenMT_gene()`, `eigenMT_batch()`)
- Data-preparation utilities — genotype PCA (`compute_geno_pcs()`), feature PCA (`compute_feature_pcs()`), rank-inverse normal transform (`rint()`)
- bigFeatures S3 class — constructor `bigFeatures()` that converts an input matrix to a disk-backed FBM
- Results are written as partitioned Apache Parquet datasets (via `arrow`)

## Repository Layout

```
R/
  bigQTL.R          # run_conditional_eqtl(), process_gene(), run_stepwise(),
                    #   run_allbutone(), input validation & sample-matching checks
  bigFeatures.R     # bigFeatures() constructor (uses bigstatsr::as_FBM())
  helpers.R         # add_snps_to_covariates(), get_cis_snps(), get_feature_indices()
  eigenMT.R         # lw_shrink_geno(), count_eigenvalues(), eigenMT convenience helpers
  prep_helpers.R    # compute_geno_pcs() (uses bigsnpr::snp_autoSVD()), compute_feature_pcs(), rint()
  bigQTL-package.R  # package-level imports / namespace hints
vignettes/
  articles/
    eigenMT_validation.Rmd  # implementation comparison of R vs Python eigenMT
tests/
  testthat.R         # testthat entrypoint (tests/ contains unit tests)
DESCRIPTION
NAMESPACE
```

## Environment Setup

### System requirements

- R ≥ 4.1
- Suggested system libraries (Debian/Ubuntu): `libhdf5-dev`, `libcurl4-openssl-dev`

### Install R dependencies

```r
# Install runtime dependencies listed in DESCRIPTION
install.packages(c("bigstatsr", "bigsnpr", "arrow"))

# Recommended for vignette/validation
install.packages(c("reticulate", "knitr", "rmarkdown"))
```

### Python environment (eigenMT comparison only)

The vignette `vignettes/articles/eigenMT_validation.Rmd` compares the R eigenMT implementation against the original Python version (an implementation comparison intended for validation and exploration, not a unit test). To reproduce the Python environment used in that vignette via `reticulate`/conda:

```r
library(reticulate)
reticulate::conda_create(
  "eigenMT",
  packages = c("numpy", "scipy", "scikit-learn", "pandas"),
  python_version = "3.8"
)
```

## Building & Checking the Package

From the repository root (shell):

```bash
R CMD build .
R CMD check bigQTL_*.tar.gz
```

From an R session (recommended for development):

```r
devtools::document()   # regenerate NAMESPACE and Rd files from roxygen blocks
devtools::build()
devtools::check()
```

Note: Documentation examples should be added as roxygen examples adjacent to the function source (in the R/ files). Rd files are generated automatically from roxygen during `devtools::document()`.

## Tests & Vignettes

- The package now contains a testthat test suite (see `tests/` and `tests/testthat.R`) and CI runs the test suite as part of checks.
- The `eigenMT_validation.Rmd` vignette performs an implementation comparison to the Python eigenMT (useful for validation and investigation). It is not a unit test and often contains interactive or long-running chunks; treat it as a reproducible comparison/validation artifact.
- CI: `.github/workflows/R-CMD-check.yaml` runs `R CMD check` on pushes and PRs and will exercise tests and example code as configured.

## Key Conventions & Notes

- File-backed matrices (FBM): genotype and feature data are intended to live on disk as `bigstatsr::FBM` objects — avoid materializing large matrices into memory.
- `bigFeatures()` constructor: converts an input matrix to an FBM using `bigstatsr::as_FBM()`.
- Genotype PCA: `compute_geno_pcs()` calls `bigsnpr::snp_autoSVD()` under the hood (handles long-range LD removal and related preprocessing).
- Sample ID matching: `design_base` row names must be meaningful sample IDs. The code matches these against `bigsnp$fam$sample.ID` and `bigfeatures$rowData$sample_name`; the package stops with an informative error if IDs are missing.
- Parallelism:
  - `ncores` — parallelism for within-gene, per-SNP operations (default 1)
  - `ncores_genes` — parallelism across genes (default 1)
- Output format: results are written as partitioned Parquet datasets. Typical layout:
  - `output_dir/stepwise/gene=<GENE>/part-*.parquet`
  - `output_dir/allbutone/gene=<GENE>/part-*.parquet`
- RINT (rank-inverse normal transform): default behavior in the `run_conditional_eqtl` pipeline is to apply RINT per feature (controlled by `do_rint`).

## Adding New Functions / Development Workflow

1. Add or edit files under `R/` (group related functions together).
2. Document public functions with roxygen2 tags in the function's source file: `@param`, `@return`, `@export`, and include examples in the roxygen `@examples` block next to the function. Use `@keywords internal` for helpers.
3. Run `devtools::document()` to update `NAMESPACE` and generate Rd files from roxygen blocks.
4. Add unit tests under `tests/testthat/` and ensure they run via `devtools::test()` locally.
5. Use `devtools::check()` locally before opening PRs.

## Contact / Issues

- Repo: https://github.com/ImmuneAxisa/bigQTL
- Please open issues for bugs, feature requests, or questions about data layout and expected inputs.
