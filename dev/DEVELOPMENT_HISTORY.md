# bigQTL — Development History

A milestone-level narrative of how the package evolved, derived from `git log`. Most
feature work was done via Copilot-agent PRs (`Initial plan` → implementation →
`Merge pull request #N` triplets); this doc collapses those into their outcomes rather
than listing every commit. Dates reflect the commit history in this repo.

## 1. Prototype: script-level QTL logic (2026-02-27)

Started as a loose collection of R functions, not yet a package:
`6a4f28c` first commit, `e5b6a77` PCA helpers, `2fc73ff`/`0143baf` fixed sample
alignment so `ind.row` is derived from `design_base` row names rather than passed in
separately (PR #1, `0a25634`) — the origin of the sample-ID-matching convention that
persists throughout the codebase today.

`81a2fbd` added `prep_helpers.R` (genotype/phenotype PCA, FBM-based, dropping the
`RNOmni` dependency in favor of an in-house `rint()`) (PR #2). `2b41fa2` moved
SNP-conditioning logic into `test_snps_with_indices()` (PR #3).

## 2. eigenMT port (2026-02-27 – 2026-03-02)

`1d1081a` added `eigenMT.R`, an R port of the Davis et al. 2016 eigenMT method (PR #4).
Validation against the original Python implementation started immediately
(`d865bda`, PR #5) but initial cross-checks failed: `34a4cf4` (2026-03-02) recorded
"Everything runs with the real eigenMT test data. But the R results are very different,
something wrong with the implementation" — the shrinkage-estimator mismatch that would
take several more passes to resolve (see §5).

## 3. Package restructuring (2026-02-28)

`ce552b6` / `c8c0d8b` turned the script collection into a proper R package (DESCRIPTION,
NAMESPACE, `R/`, roxygen docs) (PR #6). `dedf68e` moved the eigenMT validation notebook
into `vignettes/articles/` as a pkgdown-only article (not built/checked as part of
`R CMD check`) — establishing the pattern later reused for this review's new vignette
decisions. `f861492` added the first CI (`R-CMD-check.yaml`) and Copilot agent
instructions (PR #7); `326df4e` added the initial `testthat` suite. `c7aff47` added the
first comprehensive test suite for the main conditional-QTL entry point (PR #8).

## 4. Nomenclature harmonization: gene/feature → pheno (2026-03-04 – 2026-03-19)

Early code used "gene"/"feature" terminology (`bigFeatures()`, `run_conditional_eqtl()`,
`compute_feature_pcs()`, `ncores_genes`). `32c9974` harmonized this to the current
pheno-based naming (`bigPheno`, `run_conditional_qtl`, `compute_pheno_pcs`,
`ncores_phenos`) and added the `bigQTL()` convenience wrapper (PR #9/#10). The same
period added operational safeguards: `967764d` introduced `max_steps`, the `sd(y) == 0`
skip check, and `verbose` (PR #11); `177c694` added `min_snps`, `marginalQTL()`, and
moved conditioning validation into `process_pheno()` (PR #12, merged 2026-04-12).

Note: `.github/copilot-instructions.md` was written before this rename and was never
updated afterward — it still refers to `run_conditional_eqtl()`/`bigFeatures()`/
`compute_feature_pcs()`, none of which exist anymore. It was removed as part of this
review in favor of `CLAUDE.md`.

## 5. eigenMT shrinkage rewrite: OAS → Ledoit-Wolf 2004 (2026-04-14 – 2026-04-21)

This is where the "R results are very different" problem from §2 finally got resolved.
`4f21f1f` replaced the original OAS-based shrinkage with the analytical Ledoit-Wolf 2004
estimator (`lw_shrink_geno()`), matching `sklearn.covariance.LedoitWolf`'s approach of
shrinking the *raw genotype* covariance rather than a correlation matrix, and added an
optional `nlshrink`-based method for an even closer match to the Python reference
(PR #14). `0b73771` followed up addressing code review: improved the perfect-LD/
monomorphic-window edge case and switched to machine-epsilon comparisons for numerical
stability. `2728b0c` removed stale references to the old `lw_shrink_cor()` function from
the vignette (PR #16, 2026-04-21).

Also in this window: auto-documentation CI (`286b73d`/`bd30d17`/`25dff13`, PR #13,
2026-03-20/2026-03-23) — a `document` job that runs `devtools::document()` and commits
`man/`/`NAMESPACE` changes back to PR branches before `R CMD check` runs.

## 6. Stepwise conditioning fix (2026-04-23, PR #17)

The user identified a bug in an early stepwise-conditioning prototype: a conditional
regression step that failed the significance threshold was still being appended to
`stepwise_tables` before the threshold check ran, so the final results included one
"step too far" — a conditional test with nothing significant left to report.

`1e9dded` (2026-04-23) fixed `run_stepwise()` to check `new_lead$pvalue < pval_threshold`
*before* appending the step's results, and added a clarifying comment in
`run_allbutone()` documenting why its last-step reuse optimization depends on this
ordering. `8d3e217` added `set.seed(42)` for reproducibility in the new regression test
covering this exact case. Merged as PR #17 (`de9efac`, 2026-05-05).

**This bug is fixed as of the current `main`.** See `dev/CODE_REVIEW.md` for the
detailed before/after and the regression test that locks it in.

## 7. Housekeeping (2026-09-29)

`286210c` minor doc correction (`man/eigenMT_gene.Rd` wording, `.gitignore` additions
for cache/OS files); `cc64aef` added `.devcontainer/` config for a reproducible R + r2u
development environment. This review session (dev docs, code review, vignette
validation, README/CLAUDE.md rewrite) follows directly from that devcontainer setup.

## Commit reference

| Milestone | Commit(s) | Date |
|---|---|---|
| First commit | `6a4f28c` | 2026-02-27 |
| Sample-ID alignment fix | `2fc73ff`, `0a25634` | 2026-02-27 |
| eigenMT first port | `1d1081a` | 2026-02-27 |
| Package restructuring | `c8c0d8b` | 2026-02-28 |
| First CI + test suite | `f861492`, `326df4e` | 2026-02-28 |
| gene → pheno rename, `bigQTL()` wrapper | `32c9974` | 2026-03-04 |
| `max_steps`/`min_snps`/`marginalQTL()` | `967764d`, `177c694` | 2026-03-04 / 2026-03-19 |
| Auto-document CI | `286b73d` | 2026-03-20 |
| eigenMT LW2004 rewrite | `4f21f1f`, `0b73771` | 2026-04-14 |
| Stepwise "one step too far" fix | `1e9dded` (PR #17, `de9efac`) | 2026-04-23 / 2026-05-05 |
| Devcontainer | `cc64aef` | 2026-09-29 |
