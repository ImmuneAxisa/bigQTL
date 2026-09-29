# bigQTL — Code Review (2026-09-29)

Correctness-focused review of `R/*.R`. No code was changed as part of this review —
findings only. Baseline: `devtools::test()` → **69/69 passing**; `R CMD check` →
**0 errors, 0 warnings, 2 NOTEs** (both minor, see below).

**Follow-up (2026-09-29, same day): §2.1 and §2.3 fixed, and both priority gaps from
§3 closed.** `devtools::test()` → **111/111 passing** (up from 69); `covr::package_coverage()`
→ **81.10%** overall (up from 60.24%), `R/bigQTL.R` **96.42%** (up from 73.8%);
`R CMD check` → **0 errors, 0 warnings**, the `cov2cor` NOTE is gone. See the status note
at the top of each affected section below for what changed and what's still open.

## 1. Resolved: the stepwise "one step too far" bug

**Status: fixed, on `main`, with a regression test.** This was the primary thing this
review was asked to re-check.

The original prototype's `run_stepwise()` appended each conditional regression's result
table to `stepwise_tables` *before* checking whether its lead SNP passed
`pval_threshold`:

```r
# pre-fix (superseded)
stepwise_tables[[length(stepwise_tables) + 1]] <- results_step   # always appended

new_lead <- results_step[which.min(results_step$pvalue), ]
if (new_lead$pvalue < pval_threshold) {
  ...
} else {
  break   # too late — results_step is already in stepwise_tables
}
```

So the loop always computed and kept one extra conditional test beyond the last real
signal — exactly the "one step too far" behavior described.

Commit `1e9dded` (2026-04-23, PR #17) moved the append inside the `if` branch
(`R/bigQTL.R:489-494`):

```r
if (new_lead$pvalue < pval_threshold) {
  ...
  stepwise_tables[[length(stepwise_tables) + 1]] <- results_step   # only on pass
  conditioning_snps <- c(conditioning_snps, new_lead$snp)
  step <- step + 1
} else {
  ...
  break   # results_step discarded, never added
}
```

`max_steps` enforcement is also correct: the check (`bigQTL.R:464-469`) happens *before*
running the next regression, so hitting the cap never triggers an extra test either.

This is locked in by `test-run_conditional_eqtl.R:382-427`
("`run_stepwise` excludes non-passing conditioning steps from `stepwise_tables`"),
which forces a passing step-0 lead and a failing step-1 lead and asserts
`length(stepwise_tables) == 1`. **No further action needed here.**

One coupling worth being aware of when touching this code again: `run_allbutone()`'s
optimization for the last independent SNP (`bigQTL.R:539-549`) reuses
`stepwise_tables[[length(stepwise_tables)]]` on the assumption that it was computed by
conditioning on exactly `conditioning_snps[-length(conditioning_snps)]`. That assumption
only holds *because* of this fix — a regression here would silently corrupt
`run_allbutone()`'s output for the last SNP rather than error. The comment at
`bigQTL.R:541-544` documents this, which is good, but it's a fragile implicit contract
between two functions and there's no test that would catch it breaking again (see §3).

## 2. Other correctness findings

### 2.1 `eigenMT_gene`: missing `cov2cor` import (R CMD check NOTE)

**Status: fixed.** Added `@importFrom stats cov2cor` to `eigenMT_gene()`
(`R/eigenMT.R`) and regenerated `NAMESPACE`. `R CMD check` no longer NOTEs this.

`R/eigenMT.R:216` calls `cov2cor()` (used in the `"nlshrink"` branch) without importing
it from `stats`. `R CMD check` NOTEs this:

```
eigenMT_gene: no visible global function definition for 'cov2cor'
Consider adding importFrom("stats", "cov2cor") to your NAMESPACE file.
```

Not a runtime bug (`stats` is always attached), but it's a real NAMESPACE gap — if
`bigQTL` is ever used via `::` from a context where `stats` isn't attached in the search
path in the way roxygen assumes, or if this pattern is copied elsewhere, it could bite.
Fix is a one-line `@importFrom stats cov2cor` tag on `eigenMT_gene()` plus
`devtools::document()`.

### 2.2 `lw_shrink_geno()` — verified numerically against scikit-learn

Cross-checked independently during this review (not just read-through): generated a
30×12 random matrix (with one monomorphic column to exercise the degenerate-window
branch), ran `lw_shrink_geno()` and, via the `eigenMT` conda env, `sklearn.covariance.
LedoitWolf().fit()`, and compared the resulting correlation matrices. Excluding the
monomorphic row/column, **max absolute difference was 2.8×10⁻⁴** — the analytical LW2004
formula (`R/eigenMT.R:59-75`) is implemented correctly and matches the Python reference
this closely. No issue found here; the degenerate-case handling
(monomorphic window → all-ones matrix at `eigenMT.R:55-57`; near-zero-variance column →
forced `sd = 1` at `eigenMT.R:81-82`) is sound and well-commented.

### 2.3 `eigenMT_validation.Rmd` — execution results and a dead-code bug

**Status: dead code removed.** Deleted the `| lw (cvCovEst): %s` /
`comparison$m_eff_lw` / `comparison$rel_diff_lw` references from the `compare` chunk's
summary messages (the columns were never created, so these branches were always dead —
see below). The vignette's real comparison (Python vs. `"basic"` vs. `"nlshrink"`) is
unaffected; this only removed output that never printed anything in the first place.
Not re-executed against the real 373-sample eigenMT dataset as part of this fix (that
requires the conda env from the original review run); confirmed instead via
`knitr::purl()` + `parse()` that the edited chunk is still syntactically valid.

Executed end-to-end in this review (fresh conda env, cache cleared to force real
re-computation, not a cache replay) on the real 373-sample / 218,950-SNP / 15,079-gene
eigenMT example dataset, comparing Python eigenMT against both R shrinkage methods on
the first 20 genes:

```
M_eff means — Python: 5302.9 | basic: 5353.1 | nlshrink: 5302.9
Relative diff vs Python — basic: median=1% max=2%
                       nlshrink: median=0% max=0%
All 20 genes: basic within 20% of Python
```

This confirms the vignette's headline claim: `"nlshrink"` matches Python almost exactly
(0% median/max relative difference), and the dependency-free `"basic"` method tracks it
closely (1-2% relative difference) — both well inside the vignette's own 20% tolerance
gate. Rendered HTML: `vignettes/articles/eigenMT_validation.html` (gitignored build
artifact, not committed).

Not R package code, but found while validating the vignette (see task 3 of this review).
The `compare` chunk builds `comparison` with columns `m_eff_py`, `m_eff_basic`,
`m_eff_nlshrink`, `rel_diff_basic`, `rel_diff_nlshrink` (lines 377-386), but later
references `comparison$m_eff_lw` and `comparison$rel_diff_lw` (lines 398, 413), which
were never created. `data.frame$nonexistent_col` returns `NULL` rather than erroring, and
`all(is.na(NULL))` is vacuously `TRUE`, so the `if (!all(is.na(...)))` guards are always
`FALSE` — the vignette runs to completion but silently never prints the intended
"lw (cvCovEst)" comparison line. Leftover from an earlier draft that compared against a
third method. Harmless as-is (produces no wrong output, just missing output), but worth
cleaning up since it's dead/misleading code.

### 2.4 Sample-ID ordering — a load-bearing, type-unenforced invariant

**Status: now covered by a test.** `test-bigQTL.R` ("bigQTL is invariant to sample-order
permutation between design_base and bigsnp/bigpheno") constructs `design_base` with a
genuinely shuffled row order relative to `bigsnp`/`bigpheno` and asserts the resulting
stepwise output is identical (up to a handful of rows with mathematically-undefined
coefficients — see the test's comment) to the unpermuted run. The invariant itself was
not changed; this closes the "no test would catch it breaking" gap called out below.

Every entry point (`run_conditional_qtl`, `test_snps_with_indices`,
`add_snps_to_covariates`, `compute_geno_pcs`, `compute_pheno_pcs`) independently derives
row order via `match()`/`which(... %in% ...)` against `bigsnp$fam$sample.ID` /
`bigpheno$rowData$sample_name`, keyed off `rownames(design_base)`. Traced through the
full call graph in this review — it's currently **consistent everywhere**, including the
`bigQTL()` wrapper's `pheno_pcs[keep_ids, ]` / `geno_pcs[keep_ids, ]` realignment before
`cbind()`-ing PCs onto `design_base` (`bigQTL.R:196,204,206`), which works because
name-based row indexing reorders to match `keep_ids`, and `keep_ids <- rownames(design_base)`
preserves `design_base`'s existing row order. This isn't verified by any test, though
(no test constructs `bigsnp`/`bigpheno` with samples in a *different* order than
`design_base` to confirm realignment actually happens rather than accidentally lining up
because test fixtures always use matching order already — see §3). Not a bug today, but
it's exactly the kind of invariant that breaks silently (wrong sample assigned to wrong
genotype row) rather than loudly, so it deserves a dedicated test.

## 3. Test coverage gaps (measured with `covr::package_coverage()`)

**Status: both priority gaps closed (2026-09-29).** New coverage overall: **81.10%**
(was 60.24%). By file: `R/bigQTL.R` 73.8% → **96.42%**, `R/prep_helpers.R` 5.0% →
**76.67%** (picked up as a side effect of exercising `bigQTL()`), `R/bigFeatures.R`
94.4%, `R/helpers.R` 83.8% (unchanged), `R/eigenMT.R` 33.3% (unchanged — third bullet
below, not one of the two gaps prioritized here, still open).

Closing these gaps needed genuinely correlated genotype/phenotype data (existing
fixtures were all random/uncorrelated), so this pass added a reusable fixture:
`dev/generate_dummy_qtl_data.R` (seeded) simulates 5 phenotypes on one chromosome with
3/2/1/0/0 independent causal signals plus LD-friend SNPs (r² 0.2–0.9) and a background
scaffold (for realistic PC behavior), validates recovery through the package's own
`run_stepwise()`/`bigQTL()` before freezing it, and saves it as internal
`dummy_qtl_data` (`R/sysdata.rda`). `tests/testthat/helper-dummy-data.R` reconstructs
`bigsnp`/`bigpheno` from it per test (FBMs don't survive serialization, so only raw
matrices + ground truth are stored). See `dev/dummy_data_association_plots.png` for a
visual check (one line per stepwise-conditioning step, one panel per phenotype).

- ~~**`bigQTL()` — the package's namesake all-in-one wrapper — has zero test
  coverage.**~~ **Fixed:** `test-bigQTL.R` now runs `bigQTL()` end-to-end on
  `dummy_qtl_data`, asserts it recovers the correct number of independent signals per
  phenotype and identifies the true causal SNP as each step's lead, and (closing §2.4's
  gap too) asserts identical results under a genuine sample-ID permutation between
  `design_base` and `bigsnp`/`bigpheno`.
- ~~**`run_allbutone()` has zero test coverage of its actual body.**~~ **Fixed:**
  `test-run_allbutone.R` drives it directly on phenotypes with 3 and 2 real stepwise
  hits (satisfying `process_pheno`'s `length(conditioning_snps) > 1` gate for real, not
  by chance), asserts every independent SNP is recovered as a significant hit in its
  own all-but-one table, and directly checks that the fragile last-SNP-reuse
  optimization from §1 (`bigQTL.R:539-549`) matches a freshly-computed all-but-one
  result — the exact regression test that section flagged as missing.
- `eigenMT_gene()`/`eigenMT_batch()` (the FBM-facing, multi-window-splitting logic) are
  **still untested** — only their underlying pure functions (`lw_shrink_geno`,
  `count_eigenvalues`, `eigenMT_correct`) have unit tests (`test-pure-functions.R`).
  Multi-window behavior (a gene with cis-SNPs spanning more than one `window`-sized
  block) is unverified. Not one of the two gaps prioritized in this pass; `dummy_qtl_data`
  could be reused for this (e.g. lower `eigenmt_window` below a phenotype's cis-SNP
  count to force multi-window splitting) but no test was written.
- `compute_geno_pcs()`/`compute_pheno_pcs()` are now exercised indirectly via
  `test-bigQTL.R` (hence the `R/prep_helpers.R` coverage jump), but have no dedicated
  unit tests of their own (e.g. `exclude_snp_names`/`exclude_pheno_names`,
  `n_top_phenos` truncation, the `keep_ids` partial-match warning path). Still open.

## 4. R CMD check summary

Original state:

```
0 errors ✔ | 0 warnings ✔ | 2 notes ✖
❯ checking for hidden files and directories ... NOTE
  Found the following hidden files and directories: .devcontainer
  (fixed as part of this review — added to .Rbuildignore)
❯ checking R code for possible problems ... NOTE
  eigenMT_gene: no visible global function definition for 'cov2cor'
  (see §2.1 — left for the author to fix, requires an R/ code change)
```

**Status (2026-09-29): both NOTEs resolved.** `devtools::check()` now reports
**0 errors, 0 warnings, 0 notes** attributable to package code (a `.claude` hidden-dir
NOTE can appear in this sandbox — that's this AI coding session's own local tool config,
not part of the package, and won't appear in a normal checkout).
