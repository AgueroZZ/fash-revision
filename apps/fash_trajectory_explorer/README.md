# FASH trajectory explorer

A local two-tab Shiny app over all 1,009,173 time-specific-PC gene-variant pairs.

- **Trajectory explorer** — one pair at a time: observed estimates, the BF-adjusted IWP1 FASH posterior, the BF-adjusted predictive-SD mixture FASH-linear posterior, an optional Strober-style weighted trend, and four evidence cards carrying both FASH lfdr versions and both Strober statistics.
- **All-pair browser** — the whole result table: rank by any method, filter on per-method significance, keep one variant per gene, read live counts for the current view, and click a row to plot it.

## Build the compact search index

Run this command from the workflowr project root (`coderepo-local`):

```bash
Rscript --vanilla apps/fash_trajectory_explorer/build_explorer_cache.R
```

The cache builder validates key alignment across the two FASH fits, the gene map, and both Strober result tables before writing:

```text
output/revision_simulations/internal/fash_trajectory_explorer_mixture_predstep1_penalty10/explorer_index.rds
```

Both tabs read only that cache plus the two retained fits. Nothing in the app writes to `output/`, and the browser tab needs no cache rebuild.

## Launch the app

Use the launcher. It can be run from anywhere and finds the project root itself:

```bash
./apps/fash_trajectory_explorer/run_explorer.sh
```

It is idempotent: if the app is already serving it reuses that instance and just opens the browser, and if the port is free it starts the app detached, waits until the app actually answers, and then opens the browser. It never fails with "address already in use", and it refuses to touch the port if some unrelated process holds it.

```bash
./apps/fash_trajectory_explorer/run_explorer.sh status    # running or not, and on which pid
./apps/fash_trajectory_explorer/run_explorer.sh stop      # stop it (leaves foreign processes alone)
./apps/fash_trajectory_explorer/run_explorer.sh restart
./apps/fash_trajectory_explorer/run_explorer.sh log       # last 40 lines of the startup log
```

`FASH_EXPLORER_PORT` overrides the port (default 7421) and `FASH_EXPLORER_WAIT` the startup timeout (default 180 s). The startup log goes to `$TMPDIR/fash_explorer_<port>.log`; nothing is written inside the repository.

The stable local URL is <http://127.0.0.1:7421/>. The workflowr revision index links to this address, so leave `FASH_EXPLORER_PORT` unset if you want that link to work. A static HTML page cannot start an R process, so the link only resolves while the app is running — if it does not load, run the launcher.

The equivalent raw command, if you want the app in your current R session instead of detached:

```bash
Rscript --vanilla -e 'shiny::runApp("apps/fash_trajectory_explorer", host = "127.0.0.1", port = 7421, launch.browser = TRUE)'
```

The app loads the retained BF-adjusted IWP1 fit and a compact BF-adjusted FASH-linear mixture fit at startup. The compact object omits the full likelihood matrix but retains the unit-specific posterior mixture weights required for trajectory prediction. A warm-cache startup smoke check reached the listening state in about 21 seconds; memory use was not remeasured. Posterior trajectories are computed only for selected pairs and cached for the life of the app process, so the first click on a pair takes a few seconds and every later visit is instant.

Lead-variant thinning costs about 0.3 s the first time a ranking metric is used and is then memoized; a search over all 1,009,173 pairs returns in under 0.2 s.

Run the helper tests after changing anything in `explorer_helpers.R`:

```bash
Rscript --vanilla apps/fash_trajectory_explorer/test_explorer_helpers.R
```

## Tab 1: trajectory explorer

- Search accepts an HGNC symbol, Ensembl gene ID, or rsID, and updates as you type. Two terms are combined, so `GPR78 rs4583742` requires both to match. Exact matches rank first, followed by prefix and substring matches; within a match class, pairs with a smaller lfdr under either FASH model rank first. At most 200 pairs are returned.
- The ◀ and ▶ buttons step through the current result list without leaving the plot.
- Evidence cards show the headline statistic large (BF lfdr for FASH, p-value for Strober) and its companion small (raw lfdr, eFDR). The chip is the discovery call, so a statistic can never appear next to a contradicting label.
- **The weighted linear trend layer is an approximation.** Strober et al. test genotype-level data. This layer is an inverse-variance weighted least-squares line through the same time-specific-PC `beta_hat(t)` estimates that are drawn as points, is dashed, and is off by default. The exported figure repeats the caveat in its caption.
- Posterior curves are always the BF-adjusted fits, in both tabs.
- The on-screen figure omits the title and subtitle because the card header and the evidence cards already carry them; the downloaded PNG includes both so it stands alone. The per-pair CSV holds all four series on the plotting grid plus a commented evidence block.

## Tab 2: all-pair browser

- **FASH lfdr version** switches the IWP1 and FASH-linear statistic and call columns together (BF-adjusted or raw). Strober has no such distinction.
- **Rank by** sets the table order and, unless overridden, the metric used to pick each gene's lead variant. All offered statistics are smaller-is-more-significant.
- **Significance filters** are three-state per method (any / significant / not significant), so cross-method sets are expressible directly: "IWP1 significant, both Strober tests not significant" isolates FASH-only discoveries.
- **Disagreement** covers the exclusive-or cases that three-state filters cannot express: IWP1 versus FASH-linear, Strober linear versus nonlinear, and FASH (either) versus Strober (either).
- **Statistic thresholds** (collapsed) cap each statistic. A blank box means no restriction.
- **Keep one variant per gene** thins to each gene's most significant variant, breaking ties on ascending pair key. This matches `select_minimum_lfdr_variant_per_gene()` in `code/revision_simulations/internal/one_variant_per_gene_refit/`. Thinning runs on the full index *before* any other filter, so a gene's lead variant is a property of the gene over all its tested variants. With the lead metric and the significance filter on the same method, the view is exactly that method's gene-level discovery set.
- **Filtering and thinning never recalibrate FDR.** Every significance decision is the full-data cumulative-lfdr FDR-0.05 call over all 1,009,173 pairs, or the reported Strober eFDR at 0.05. The counts in the summary tiles describe the current view under those fixed calls; they are not FDR-controlled results for the subset. Refitting on a thinned set is a separate analysis, and that is what `one_variant_per_gene_refit` does.
- The **Calls** column packs the four discovery calls into four dots (IWP1 FASH, FASH-linear, Strober linear, Strober nonlinear, left to right; filled means significant). It is rendered from a 4-bit code in the browser so a million-row view never has to build a million HTML strings.
- The table has no DataTables search box on purpose: it would filter the page without updating the summary tiles. Use the sidebar text filter, which does both.
- **Download current view (CSV)** writes the full 18-column index rows for exactly the pairs in view, in view order, plus the call code and the ranking column name.
