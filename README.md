# bigQTL

bigQTL is an R package for scalable conditional cis-QTL analysis that is built to work with file-backed matrices (FBM) from the bigstatsr / bigsnpr ecosystem. It provides fast within-gene and across-gene workflows for stepwise conditional QTL mapping, all-but-one conditioning, and eigenMT multiple-testing correction — all while keeping large genotype and phenotype matrices on disk.

Key goals:
- Scale cis-QTL discovery to large sample sizes using disk-backed FBMs
- Provide reproducible, partitioned Parquet outputs for downstream analysis
- Offer eigenMT-based multiple-testing correction with a Ledoit–Wolf shrinkage estimator

Repository: https://github.com/ImmuneAxisa/bigQTL

Table of contents
- Features
- System requirements
- Installation
- Quickstart examples
- Output format
- Conventions & best practices
- Vignettes & validation
- Development & contributing
- License & citation
- Contact

Features
- Conditional QTL analysis
  - bigQTL(): all-in-one wrapper with automatic PC computation
  - run_conditional_qtl(): stepwise conditioning and all-but-one conditioning
  - process_pheno(), run_stepwise(), run_allbutone()
- eigenMT multiple-testing correction
  - lw_shrink_geno(), count_eigenvalues(), eigenMT_gene(), eigenMT_batch()
- Data-preparation utilities
  - compute_geno_pcs(), compute_pheno_pcs(), rint()
- Helper utilities
  - get_pheno_indices(), get_snp_indices(), get_cis_snps(), add_snps_to_covariates()
- bigPheno S3 class for wrapping phenotype/expression matrices in an FBM
- Output written as partitioned Apache Parquet datasets (via arrow)
- Parallelism controls:
  - `ncores` — parallelism within-phenotype (SNP-level)
  - `ncores_phenos` — parallelism across phenotypes

System requirements
- R >= 4.1
- Suggested system libraries for package installation:
  - libhdf5-dev
  - libcurl4-openssl-dev

R package dependencies (declared in DESCRIPTION)
- bigstatsr
- bigsnpr
- arrow
(Development / vignette only)
- reticulate, knitr, rmarkdown

# Installation

Install system libraries (example for Debian/Ubuntu):
```bash
sudo apt-get update
sudo apt-get install -y libhdf5-dev libcurl4-openssl-dev
```

Install R dependencies:
```r
install.packages(c("bigstatsr", "bigsnpr", "arrow"))
# For vignette / validation:
install.packages(c("reticulate", "knitr", "rmarkdown"))
```

Install bigQTL from GitHub:
```r
# one of:
# with remotes
remotes::install_github("ImmuneAxisa/bigQTL")

# or clone and install locally
git clone https://github.com/ImmuneAxisa/bigQTL.git
cd bigQTL
R CMD build .
R CMD INSTALL bigQTL_*.tar.gz
```

Python environment for eigenMT validation (optional)
The repository includes an eigenMT validation vignette that compares the R implementation to the original Python version. To create the conda environment used in the vignette:
```r
library(reticulate)
reticulate::conda_create(
  "eigenMT",
  packages = c("numpy", "scipy", "scikit-learn", "pandas"),
  python_version = "3.8"
)
```

# Quickstart

Below are short examples demonstrating common workflows. These assume you have:
- A bigSNP object (from bigsnpr) or genotype FBM
- A bigPheno object or expression FBM
- A design matrix with sample IDs that match the genotype and phenotype objects

1) Construct a bigPheno object (example)
```r
library(bigQTL)
# Suppose expr_mat is a matrix with rows = samples, cols = phenotypes
bp <- bigPheno(expr_mat, rowData = my_rowdata)
```

2) Compute genotype PCs (recommended)
```r
pcs <- compute_geno_pcs(bigsnp, n_pcs = 10, ncores = 4)
```

3) Run all-in-one QTL analysis (computes PCs automatically)
```r
res <- bigQTL(
  bigpheno    = bp,
  bigsnp      = bigsnp,
  pheno_coord = pheno_coord_df,  # data.frame with pheno_name, chromosome, start, end
  design_base = design_matrix,   # rownames must be sample IDs
  n_pheno_pcs = 5,
  n_geno_pcs  = 5,
  output_dir  = "output/"
)
```

4) Run stepwise conditional cis-QTL analysis directly
```r
res <- run_conditional_qtl(
  bigpheno    = bp,
  bigsnp      = bigsnp,
  pheno_coord = pheno_coord_df,
  design_base = design_matrix,
  output_dir  = "output/",
  ncores      = 4,
  ncores_phenos = 8
)
# Writes partitioned Parquet files under output_dir/stepwise/pheno=<name>/
```

5) Apply eigenMT correction per phenotype or in batch
```r
# Run eigenMT over a batch of phenotypes
eigenMT_batch(bigsnp, pheno_coord = pheno_coord_df, ind.row = ind_row, ncores = 4)
```

# Notes

Output format
- Results are written as partitioned Apache Parquet datasets using the arrow package.
- Output layout:
  - output_dir/stepwise/pheno=<PHENO>/part-0.parquet
  - output_dir/allbutone/pheno=<PHENO>/part-0.parquet
- Parquet schema is column-oriented and partitioned by phenotype for easy downstream aggregation.
- Each row corresponds to a SNP test (or conditional step) with metadata: phenotype, SNP id, position, effect size, SE, p-value, q-value, conditioned_snps, step_index, etc.

Conventions & best practices
- File-backed matrices (FBM): All genotype and phenotype data should live on disk as bigstatsr::FBM objects. Avoid materializing large matrices into memory.
- Sample ID matching: design_base rownames must be meaningful sample IDs that appear in both bigsnp$fam$sample.ID and bigpheno$rowData$sample_name.
- Parallelism:
  - ncores: parallelism for inner loops (SNP-level operations)
  - ncores_phenos: parallelism across phenotypes (useful for multi-core servers)
- RINT (rank inverse normal transform): default do_rint = TRUE. Keeps phenotypes well-behaved across phenotypes.
- eigenMT: uses Ledoit–Wolf shrinkage estimator for stable correlation estimates when sample size is limited compared to number of SNPs.

Vignettes & validation
- The repository contains a vignette that validates the R eigenMT implementation against the original Python implementation:
  - vignettes/articles/eigenMT_validation.Rmd
  - Link: https://github.com/ImmuneAxisa/bigQTL/blob/main/vignettes/articles/eigenMT_validation.Rmd
- Note: Interactive chunks in the vignette are set with eval = FALSE and are intended to be run locally.


If you plan to open a pull request:
- Create a branch from main, implement your changes, add tests (where possible), update documentation and vignettes, and open a PR describing the change and rationale.

License & citation
- See the LICENSE file in the repository for licensing details.
- When using results from bigQTL in publications, please cite the repository and include appropriate method description (e.g., eigenMT and Ledoit–Wolf shrinkage references).

Acknowledgements & References
- bigstatsr / bigsnpr — file-backed matrices and genotype utilities
- eigenMT — multiple-testing correction approach for cis-QTL mapping
- Ledoit & Wolf — shrinkage covariance estimator

Contact
- Repository: https://github.com/ImmuneAxisa/bigQTL
- Issues and feature requests: file an issue in the repository
- For questions, open a discussion or issue detailing your use case and data setup

Notes
- bigQTL is designed for large-scale, production-style QTL analyses where disk-backed storage and memory efficiency are critical. If you are exploring small proofs-of-concept, in-memory approaches may be simpler — but for realistic genotype / phenotype sizes, use FBMs and the provided workflows.
