# bigQTL — Architecture Overview

Internal developer reference for how the package fits together. For the public-facing
summary see `README.Rmd`; for how it got this way see `dev/DEVELOPMENT_HISTORY.md`; for
known issues see `dev/CODE_REVIEW.md`.

## Data model

bigQTL operates entirely on **file-backed matrices (FBM)** from the bigstatsr/bigsnpr
ecosystem so that genotype/phenotype matrices larger than RAM can be analyzed without
materializing them in memory.

Two parallel container objects:

| Object | Constructor | Backing matrix | Row metadata | Column metadata |
|---|---|---|---|---|
| `bigSNP` | (from bigsnpr, not bigQTL) | `$genotypes` FBM (samples × SNPs) | `$fam$sample.ID` | `$map$marker.ID`, `$map$chromosome`, `$map$physical.pos` |
| `bigPheno` | `bigPheno()` (`R/bigFeatures.R:16`) | `$pheno` FBM (samples × phenotypes) | `$rowData$sample_name` | `$colData$pheno_name` |

`bigPheno` deliberately mirrors bigsnpr's `bigSNP` shape so the rest of the package can
treat genotypes and phenotypes symmetrically (match samples by name, look up columns by
name via `get_snp_indices()` / `get_pheno_indices()`).

**Load-bearing invariant:** every entry point re-derives row order from
`rownames(design_base)` and uses `match()` against `bigsnp$fam$sample.ID` /
`bigpheno$rowData$sample_name` to build `ind.row.snp` / `ind.row.pheno`
(`run_conditional_qtl()`, `bigQTL.R:59-73`). All downstream indexing (`y`, covariates,
conditioning-SNP columns) assumes this same sample order. Nothing enforces this at the
type level — if you add a new code path that touches `y`/`design_base`/genotype columns,
it must go through the same `ind.row.*` derivation or samples will silently misalign.

## File-by-file catalog

### `R/bigQTL.R` — orchestration + stepwise/all-but-one conditioning

- `run_conditional_qtl()` (32) — main engine. Validates inputs (sample-ID matching,
  `pheno_coord` columns, `do_allbutone` requires `do_conditioning`), fans out over
  phenotypes via `lapply`/`parallel::mclapply` (`ncores_phenos`), writes one Parquet
  partition per phenotype, then opens the result as an Arrow dataset.
- `bigQTL()` (183) — convenience wrapper: computes phenotype PCs (`compute_pheno_pcs()`)
  and genotype PCs (`compute_geno_pcs()`), appends them to `design_base`, calls
  `run_conditional_qtl()`.
- `marginalQTL()` (244) — wrapper that disables both conditioning modes
  (`do_conditioning = do_allbutone = FALSE`) for a marginal-only scan.
- `process_pheno()` (294, internal) — per-phenotype pipeline: RINT transform → cis-SNP
  lookup (`get_cis_snps()`) → marginal test (step 0) → `run_stepwise()` →
  `run_allbutone()` → write Parquet. Skips a phenotype (with a warning) if `sd(y) == 0`
  or fewer than `min_snps` cis-SNPs are found.
- `run_stepwise()` (440, internal) — greedy forward selection. Loop: test all cis-SNPs
  conditioned on the SNPs accepted so far; if the new lead SNP's p-value passes
  `pval_threshold`, keep the step and add the SNP to the conditioning set; otherwise
  discard the step and stop. Bounded by `max_steps` (checked before running the next
  regression, so hitting the cap never runs one extra test). See
  `dev/CODE_REVIEW.md` for the history of a bug here (now fixed) where a failing step
  was appended to results before the significance check.
- `run_allbutone()` (526, internal) — for each SNP in `conditioning_snps`, re-test all
  cis-SNPs conditioned on every *other* accepted SNP. The last SNP's all-but-one result
  is provably identical to the last accepted stepwise table (conditioning on
  `conditioning_snps[-length]`), so it's reused instead of recomputed.
- `test_snps_with_indices()` (593, exported) — the actual regression: optionally expands
  `design_base` with conditioning-SNP genotype columns (`add_snps_to_covariates()`), then
  calls `bigstatsr::big_univLinReg()` for a vectorized per-SNP univariate linear
  regression, adding BH-FDR via `stats::p.adjust()`.

### `R/eigenMT.R` — eigenMT multiple-testing correction

Implements Davis et al. 2016 (eigenMT): for each phenotype's cis-SNP window(s), estimate
the effective number of independent tests (`M_eff`) from the eigenvalue spectrum of the
LD correlation matrix, rather than assuming every SNP is an independent test.

- `lw_shrink_geno()` (39, exported) — analytical Ledoit-Wolf 2004 shrinkage applied to
  the *raw genotype* covariance (mirrors `sklearn.covariance.LedoitWolf`), converted to a
  correlation matrix. No package dependency beyond base R. Handles degenerate windows
  (monomorphic SNPs → all-ones matrix; near-zero-variance columns → forced sd = 1,
  independent-test fallback).
- `count_eigenvalues()` (104, exported) — how many top eigenvalues are needed to explain
  `var_thresh` (default 0.99) of total variance.
- `eigenMT_correct()` (136, exported) — `pmin(pvalue * m_eff, 1)`.
- `eigenMT_gene()` (176, exported) — splits a phenotype's cis-SNPs into disjoint windows
  (default 200 SNPs), shrinks + eigendecomposes each window, sums `count_eigenvalues()`
  across windows. `shrinkage_method` picks `"basic"` (`lw_shrink_geno()`, default) or
  `"nlshrink"` (`nlshrink::linshrink_cov()`, optional dependency, matches the Python
  reference implementation more closely — see `dev/CODE_REVIEW.md` /
  `vignettes/articles/eigenMT_validation.Rmd`).
- `eigenMT_batch()` (262, exported) — runs `eigenMT_gene()` for every row of
  `pheno_coord`, parallelized across phenotypes via `parallel::mclapply`.

### `R/helpers.R` — name/index lookups and covariate augmentation

- `add_snps_to_covariates()` (20) — appends genotype columns for named SNPs to a
  covariate data frame (used by `test_snps_with_indices()` for conditioning).
- `get_cis_snps()` (69) — filters `bigsnp$map` to SNPs on the same chromosome within
  `[start - cis_window, end + cis_window]`.
- `get_pheno_indices()` (105) / `get_snp_indices()` (131) — name → FBM column index,
  erroring on any unmatched name.

### `R/prep_helpers.R` — data preparation

- `compute_geno_pcs()` (28) — `bigsnpr::snp_autoSVD()` (handles long-range LD/outlier
  removal automatically), returns a named PC matrix.
- `compute_pheno_pcs()` (97) — operates on the phenotype FBM directly:
  `big_colstats()` to rank phenotypes by variance (avoids PCA on every column),
  `big_randomSVD()` on the top `n_top_phenos` (default 5000).
- `rint()` (162) — Blom rank-inverse-normal transform, a dependency-free replacement for
  `RNOmni::RankNorm()`.

### `R/bigFeatures.R` — `bigPheno` constructor

Validates the input matrix/rowData/colData shapes, wraps the matrix in an FBM via
`bigstatsr::as_FBM()`, and returns the S3 list described in "Data model" above.

## Pipeline / data flow

```
Prep                         QTL mapping                          MT correction
----                         -----------                          -------------
bigPheno(matrix)      \
bigSNP (from bigsnpr)  |--> run_conditional_qtl() / bigQTL() -->  eigenMT_batch()
design_base (covars)  /      (per phenotype, parallel over        (per phenotype,
pheno_coord                   phenotypes via ncores_phenos):        parallel via ncores)
                                process_pheno():                  -> data.frame(pheno_name,
                                  RINT(y) if do_rint                 n_cis_snps, m_eff)
                                  get_cis_snps()
                                  test_snps_with_indices()  (step 0)   eigenMT_correct(pvalue, m_eff)
                                  run_stepwise()             (steps 1..N)
                                  run_allbutone()            (leave-one-out)
                                  write_parquet()
                              -> arrow::open_dataset() per output type
```

Output layout: `output_dir/{stepwise,allbutone}/pheno=<name>/part-0.parquet`, opened as
Arrow datasets so results across all phenotypes can be queried lazily
(`arrow::open_dataset() |> dplyr::filter(...) |> dplyr::collect()`).

## Parallelism

Two independent knobs, safe to combine:

- `ncores` — within a single phenotype's regression, passed straight to
  `bigstatsr::big_univLinReg(..., ncores = ncores)` (and to `bigsnpr::snp_autoSVD()` /
  `bigstatsr::big_randomSVD()` in the prep helpers).
- `ncores_phenos` — across phenotypes, via `parallel::mclapply()` in
  `run_conditional_qtl()` and `eigenMT_batch()` (there the parameter is just `ncores`,
  since there's no separate "within" stage).

`mclapply`-based parallelism is fork-based (POSIX only; no Windows support).
