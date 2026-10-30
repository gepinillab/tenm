# Code review of `tenm` 0.5.1 (`R/`)

**Scope:** every file in `R/` at commit `37b4cdf`. No package code was changed. This file is the only change in this branch.

**How the findings were checked:** I read every function, then ran the workflow from the README and vignette (`sp_temporal_data` → `clean_dup_by_date` → `ex_by_date` → `bg_by_date` → `tenm_selection` → `predict`) on the bundled `abronia` data and `inst/extdata/bio` layers. Each finding marked **Confirmed** below has a minimal reproducible example (MRE) that I ran; the output shown is real output.

**Environment:** R 4.3.3 on Ubuntu 24.04 with GDAL/GEOS/PROJ. I ran everything first with the Ubuntu-packaged dependencies (terra 1.7-65, sf 1.0-15, purrr 1.0.2, tidyr 1.3.1, furrr 0.3.1, future 1.33.1, MASS 7.3-60). I then re-ran the full pipeline and the repros for B1, B2, B3, B4, B7, B9, B10, B12 and R7 against current CRAN releases (terra 1.9-50, sf 1.1-3, purrr 1.2.2, tidyr 1.3.2, dplyr 1.2.1, furrr 0.4.0, future 1.80.0). Results were identical on both stacks.

All MREs start from this shared setup, which caches the README objects:

```r
library(tenm); data("abronia"); set.seed(1)
d <- system.file("extdata/bio", package = "tenm")
tenm_mask <- terra::rast(file.path(d, "1939/bio_01.tif"))
abt  <- sp_temporal_data(abronia, "decimalLongitude", "decimalLatitude", "year", "y", "y", d)
abtc <- clean_dup_by_date(abt, by_mask = TRUE, threshold = terra::res(tenm_mask)[1],
                          raster_mask = tenm_mask, n_ngbs = 0)
abex <- ex_by_date(abtc, 0.7)
abbg <- bg_by_date(abex, buffer_ngbs = 10, n_bg = 10000)
```

Severity scale:
- **bug:** the code gives wrong results or fails on valid input.
- **risk:** the code works on the bundled example but is fragile or silently wrong under realistic conditions.
- **improvement:** performance, maintainability, or CRAN hygiene.

---

## Summary

| # | Sev. | Where | One-line summary |
|---|------|-------|------------------|
| B1 | bug | `ellipsoid_selection.R:203,240` | `proc_iter` passed by position binds to `mve`, so `mve = FALSE` and `NoOfIteration` are silently ignored |
| B2 | bug | `00_methods.R:190-196` | `predict()` stacks layers in file order, not model order, so the projection is wrong whenever the variable order differs |
| B3 | bug | `ellipsoid_omr.R:188` | Binomial p-value is off by one (`P(X > k)` instead of `P(X ≥ k)`), making it anti-conservative (0 whenever all test points are predicted) |
| B4 | bug | `bg_by_date.R:59-115` | Background is aligned to years by row position, so it can be drawn around the wrong year's occurrences and get the wrong per-year sample size |
| B5 | bug | `ex_by_date.R:87-103` | Train/test labels are assigned by row position after `split()`, so they land on the wrong rows when folder order differs from date order |
| B6 | bug | `ex_by_date.R:46-54`, `bg_by_date.R:123,142` | Cell IDs from the first year's grid are reused for every year, so background values come from wrong pixels when grids differ |
| B7 | bug | `clean_dup_by_date.R:89-94` | Fails under any `multisession` plan when an argument is a variable (`globals=` whitelist) |
| B8 | bug | `clean_dup.R:69`, `README.Rmd:242` | `by_mask = TRUE` with a non-`SpatRaster` mask silently falls back to distance thinning; the README itself does this |
| B9 | bug | `clean_dup.R:150-162` | Distance thinning removes records whose only neighbour was itself removed (over-thins chains) |
| B10 | bug | `clean_dup.R:76-77,141` | Off-mask points: all but one are dropped (`n_ngbs = 0`) or all are dropped (`n_ngbs > 0`), silently |
| B11 | bug | `bg_by_date.R:102` | `sample(x, n)` with `length(x) == 1` samples from `1:x`, giving a random wrong cell and then a crash |
| B12 | bug | `PartialROC.R:131,193` | `pROC()` ignores `rseed` and always calls `set.seed(1..n_iter)`, clobbering the user's RNG stream |
| B13 | bug | `ellipsoid_omr.R:182` | `exists("in_etest")` also finds objects in the global environment, so the call crashes when the user has an `in_etest` |
| B14 | bug | `tenm_selection.R:109-112` | `RandomPercent`, `sub_sample`, `sub_sample_size` are documented but never used |
| B15 | bug | `00_methods.R:133,151,177,203,193` | `predict()`: `plot = FALSE` errors, the bad-variable message prints `NULL`, and the default `layers_ext = NULL` errors |
| B16 | bug | `sp_temporal_data.R:60,21-23` | Tibble input crashes; `"mdy"` is documented but rejected |
| B17 | bug | `ex_by_date.R:84` | Exact duplicate records produce list-columns that break every later step |
| B18 | bug | `cov_center.R:53` | Default `vars = NULL` selects zero columns, giving a misleading "collinear" error |
| B19 | bug | `ellipsoid_omr.R:124` | Tests `fails_train_ids` instead of `fails_test_ids`, so `non_pred_test_ids` is `""` instead of `NA` |
| B20 | bug | `clean_dup.R:65,85-87` | `NA` latitude crashes (only longitude is checked); coordinate-only input crashes with `by_mask = TRUE` |
| B21 | bug | `PartialROC.R:217-234` | `proc_precision()` errors when the test data have fewer than 2 unique values, so the pROC result becomes `NA` |
| R1 | risk | `sp_temporal_data.R:69-75,114-118` | Occurrences with no layer for their year are silently matched to the nearest year (1890 → 1939, 2050 → 2016); unparseable dates are silently dropped |
| R2 | risk | `00_methods.R:143`, `tenm_selection.R` | `predict()` refits the model with its own `mve`/`level` defaults (and new MVE randomness) instead of the selected model's |
| R3 | risk | several | No validation of CRS, extent, resolution, or layer completeness across years |
| R4 | risk | `ellipsoid_selection.R:251-254,284`, `ellipsoid_omr.R:69-77` | Missing test omission is set to 0; failed models vanish silently; the output schema changes when no model passes |
| R5 | risk | `ellipsoid_omr.R:163` | Background prevalence also counts training presences, which inflates the prevalence used by the binomial test |
| R6 | risk | `ex_by_date.R:77,111`, `tdf2swd.R:41-48`, `clean_dup.R:85` | Columns are addressed by position, and `clean_dup` reorders them |
| R7 | risk | `ex_by_date.R:73`, `bg_by_date.R:130`, `PartialROC.R:63-66`, `metaras.R:13` | `eval(parse())` on terra's internal C++ module fields |
| R8 | risk | `ellipsoid_selection.R:131-177` | Parallel path: overrides the user's `future` plan, leaks options, ignores RNG misuse, so results are not reproducible |
| P1 | improvement | `bg_by_date.R:146-159` | Per-row `map_dfc` unpacking takes 85% of runtime (~9 ms per background point) |
| P2 | improvement | `ellipsoid_omr.R`, `PartialROC.R` | pROC runs for every candidate model and takes 83% of selection time |
| P3 | improvement | `clean_dup.R:148-150` | Dense n×n intersection matrix: memory grows quadratically (1.2 GB at n = 8,000) |
| P4 | improvement | `ellipsoid_projection.R:66-96` | Loads whole rasters into memory and builds a data.frame column by column |
| P5 | improvement | `ex_by_date.R:56-81`, `bg_by_date.R:93,118-137` | One raster open plus a full-table filter per file; the `buffer_distance` branch opens one copy of the layer per occurrence |
| D1-D6 | improvement | DESCRIPTION, tests, examples | Dependency and CRAN-check items (see section D) |
| M1-M7 | improvement | several | Duplicated code, very long functions, API inconsistencies, tests that check only classes |

**Deprecated dependencies:** the package no longer uses `raster`, `sp`, `rgdal`, `rgeos` or `maptools` (the only mention is a commented-out `sp::` line in `clean_dup.R:68`). The remaining dependency issues are listed in section D.

---

## A. Confirmed bugs

### B1 (bug): `ellipsoid_selection()` silently ignores `mve` and the number of pROC iterations

**Where:** `R/ellipsoid_selection.R:203` (parallel branch) and `:240` (serial branch).

**What's wrong:** both branches call

```r
tenm::ellipsoid_omr(env_data = ..., env_test = ..., env_bg = ..., cf_level = level,
                    proc = proc, proc_iter, rseed = rseed)
```

`proc_iter` is unnamed, so R binds it to the first unmatched formal argument of `ellipsoid_omr`, which is `mve`, not `proc_iter`. As a result:

- `mve` becomes the iteration count (for example 1000). That is truthy, so an MVE is always fitted, and `ellipsoid_selection(mve = FALSE)` / `tenm_selection(mve = FALSE)` never fit the classical covariance ellipsoid.
- `proc_iter` falls back to its default of 100, so `tenm_selection(NoOfIteration = 1000)` runs 100 bootstrap iterations.

This affects every calibration run and every reported p-value and AUC ratio.

**Evidence:**

```r
match.call(tenm::ellipsoid_omr, quote(tenm::ellipsoid_omr(env_data = env_data, env_test = env_test,
           env_bg = env_bg, cf_level = level, proc = proc, proc_iter, rseed = rseed)))
#> tenm::ellipsoid_omr(env_data = env_data, env_test = env_test, env_bg = env_bg,
#>     cf_level = level, mve = proc_iter, proc = proc, rseed = rseed)

trace("ellipsoid_omr", quote(cat("mve =", mve, "| proc_iter =", proc_iter, "\n")),
      where = asNamespace("tenm"), print = FALSE)
ellipsoid_selection(etr, ete, c("bio_01","bio_12"), 2, 0.975, mve = FALSE,
                    env_bg = abbg$env_bg[, c("bio_01","bio_12")], omr_criteria = .1,
                    proc = TRUE, proc_iter = 7)
#>   ellipsoid_omr got mve = 7 | proc_iter = 100

# Effect on results (same data, same seed):
#> bg_prevalence  direct mve=FALSE: 0.5212  direct mve=TRUE: 0.5004  ellipsoid_selection(mve=FALSE): 0.5004
```

**Fix:** name the argument (`proc_iter = proc_iter`) and pass `mve = mve` explicitly in both branches. Add a regression test that `ellipsoid_selection(mve = FALSE)` matches `ellipsoid_omr(mve = FALSE)`. Fixing this makes pROC about 10× slower at the documented 1000 iterations (see P2).

---

### B2 (bug): `predict()` mismatches variables and layers when the variable order differs from the file order

**Where:** `R/00_methods.R:190-196` (`layers_path` branch). The `SpatRaster` and `list` branches (`:147`, `:174`) have the same problem: they only check the layer count (`:162`), never the names.

**What's wrong:** the centroid and covariance come from `cov_center(env_data, vars = mod_vars)`, so they follow the order of `mod_vars`. The layers come from `list.files(...)[var_ids]`, so they follow alphabetical file order. `ellipsoid_projection()` then calls `stats::mahalanobis()`, which matches columns by position. Whenever `model_variables` (or a `fitted_vars` entry built from a user-ordered `vars2fit`) is not alphabetical, every pixel is scored against swapped axes. No error is raised and the map is simply wrong.

**Evidence:** the same two-variable model, listed in two orders:

```r
p1 <- predict(mod_sel, model_variables = c("bio_01","bio_12"), layers_path = dir70, layers_ext = ".tif$")
p2 <- predict(mod_sel, model_variables = c("bio_12","bio_01"), layers_path = dir70, layers_ext = ".tif$")
max(abs(values(p1) - values(p2)), na.rm = TRUE)  #> 0.9976
global(c(p1, p2), "mean", na.rm = TRUE)          #> 1.94e-01  vs  3.96e-65
# Reference: layers explicitly ordered by name, then ellipsoid_projection():
#   max diff vs p1 = 1.1e-16 (p1 is correct), max diff vs p2 = 0.9976 (p2 is wrong)
```

**Fix:** select layers by name, not by position: `slayers <- terra::rast(files)[[mod_vars]]` after setting `names()` from the file stems. In `ellipsoid_projection()`, check that `all(names(centroid) %in% names(envlayers))` and reorder with `envlayers[[names(centroid)]]`. Fail with a clear message listing any missing variables; today a missing variable gives `non-conformable arguments` from `x %*% cov`.

---

### B3 (bug): the binomial-test p-value is off by one (anti-conservative)

**Where:** `R/ellipsoid_omr.R:188-190`

**What's wrong:** `p_bin <- 1 - pbinom(test_succs, size, prob)` computes `P(X > k)`. The one-sided p-value for observing at least `k` successes is `P(X ≥ k) = 1 - pbinom(k - 1, ...)`. The current formula underestimates the p-value, and it is exactly 0 whenever every test point is inside the ellipsoid, whatever the prevalence. That is the case most likely to be selected.

**Evidence:** 5 of 5 test points inside, with a background prevalence of 0.925:

```r
#> tenm pval_bin = 1 - pbinom(k, n, p)     = 0
#> correct        1 - pbinom(k - 1, n, p) = 0.678   (binom.test(..., alternative = "greater") = 0.678)
```

**Fix:** use `stats::pbinom(test_succs - 1, size, prob, lower.tail = FALSE)`, or `binom.test(test_succs, n, prevBG, alternative = "greater")$p.value`.

---

### B4 (bug): `bg_by_date()` aligns background to years by row position

**Where:** `R/bg_by_date.R`

- `:61-70`: `paths_layers` is built per `split(tdf, tdf$layers_path)` group (sorted order), then assigned to `tdf$layer_path` row by row (original order).
- `:107, :115`: `names(cells_to_samp) <- unique(tdf$layers_path)` uses order of appearance, but `cells_to_samp` is in `split()` (sorted) order.
- `:59` vs `:84`: `df_samps` is sorted by `layer_dates`, `ddL` is sorted by `layer_path` string, and `df_samps$n_samples[z]` is indexed by position.

**What's wrong:** the result is only correct if (a) `temporal_df` is already sorted by `layers_path` and (b) the alphabetical order of folder names equals their date order. (a) holds for unmodified `ex_by_date()` output, but breaks as soon as a user filters, reorders, or binds `temporal_df`. (b) breaks for `"my"`/`"dmy"` folders such as `01-1940` and `12-1939`. When either fails, background cells sampled around one year's occurrences are extracted from another year's layers, or a year gets another year's sample size.

**Evidence 1:** the same `abex`, with rows reversed. I measured the share of each year's background that lies within `buffer_ngbs` of that year's occurrences:

```r
ok  <- bg_by_date(abex, buffer_ngbs = 2, n_bg = 1000)
ex2 <- abex; ex2$temporal_df <- ex2$temporal_df[nrow(ex2$temporal_df):1, ]
bad <- bg_by_date(ex2, buffer_ngbs = 2, n_bg = 1000)
#> sorted input: mean share 1.000 | rows reversed: mean share 0.492 (26 of 32 years affected)
```

**Evidence 2:** folders `12-1939` (1 occurrence) and `01-1940` (3 occurrences), `n_bg = 40`:

```r
#> expected 10 for 12-1939 and 30 for 01-1940; got:
#> 01-1940 12-1939
#>      10      30
```

**Fix:** never rely on positional alignment after `split()`. Carry the key along: for example, build one lookup `data.frame(layers_path, layer_dates, n_samples, template)` with `dplyr::summarise(.by = layers_path)`, iterate over it, and name results by `layers_path`. For example, `cells_to_samp <- setNames(purrr::map(keys$layers_path, ...), keys$layers_path)`.

---

### B5 (bug): `ex_by_date()` assigns train/test labels to the wrong rows

**Where:** `R/ex_by_date.R:87-103`

**What's wrong:** labels are generated per `split(years_env, years_env$layer_dates)` group, which is sorted by date. They are then `unlist()`ed and assigned positionally to `years_env`, whose rows are in file-path order from `list.files`. When those orders differ, labels land on other years' rows. The per-year guarantees in the code (a single record is always Train; two records give one Train and one Test; `train_prop` per year) are then lost.

**Evidence:** folders `12-1939` and `01-1940` (`"my"` format), with 1 and 2 occurrences:

```r
#>     date layer_dates trian_test
#> 1 01-1940  1940-01-01      Train
#> 2 01-1940  1940-01-01      Train     # the 2-record year has no Test record
#> 3 12-1939  1939-12-01       Test     # the 1-record year was meant to be Train
```

**Fix:** assign within groups, for example `years_env <- dplyr::mutate(years_env, trian_test = assign_split(dplyr::n(), train_prop), .by = layers_path)`, or use `split<-` / `unsplit()` with the same factor.

---

### B6 (bug): cell IDs from the first year's raster are used for every year

**Where:** `R/ex_by_date.R:46-54` (`cell_ids_year` computed from `tdf$layers_path[1]` only). These IDs are then used in `bg_by_date.R:78-86` (buffers), `:123` (`env_layers[cellids]`) and `:142` (`xyFromCell(r1, ...)`).

**What's wrong:** the IDs only identify the same location if every year's rasters share the same extent and resolution. Nothing checks this. If one year's grid differs (a different source, a crop, or a different resolution), background values come from different pixels than the coordinates reported next to them. Occurrence values are not affected, because `ex_by_date` extracts them by coordinates.

**Evidence:** I cropped the 1940 layers by 3 columns (same CRS and resolution) and ran the pipeline on 1939 and 1940:

```r
#> 1940 background: reported bio_01 vs value at the reported coordinates
#>   reported at_xy
#> 1      133   196
#> 2      165   128
#> 3      147   141
#> ...
#> mismatching rows: 14 of 14
```

**Fix:** validate the geometry once, in `sp_temporal_data()` or `ex_by_date()`: `terra::compareGeom()` across one layer per year folder, with an informative error. Alternatively, compute cell IDs per year from that year's template and extract background by coordinates.

---

### B7 (bug): `clean_dup_by_date()` fails under `future::plan(multisession)`

**Where:** `R/clean_dup_by_date.R:89-94`

**What's wrong:** `furrr_options(globals = c("this_species", "df_occs_date", "df_occs_dateL"))` turns off automatic global detection. `threshold`, `raster_mask`, `by_mask` and the other arguments are lazy promises that point at the caller's environment, which is not exported. Any non-literal argument therefore fails on the workers. Even if forced, a `SpatRaster` cannot be serialized to a worker without `terra::wrap()`. The README's recommended pattern of setting `plan("multisession")` early in a script makes this function fail.

**Evidence:**

```r
future::plan("multisession", workers = 2)
th <- 0.1
clean_dup_by_date(abt, threshold = th)
#> Caused by error: ! object 'th' not found
clean_dup_by_date(abt, threshold = 1, by_mask = TRUE, raster_mask = tenm_mask)
#> Caused by error: ! object 'tenm_mask' not found
```

**Fix:** drop the `globals =` whitelist and `force()` the arguments, or pass them via `.f`'s arguments. Pass `terra::wrap(raster_mask)` and `terra::unwrap()` it inside the worker. The work per group is small, so `purrr::map()` is probably enough here.

---

### B8 (bug): `by_mask = TRUE` with an invalid `raster_mask` silently uses distance thinning, and the README does this

**Where:** `R/clean_dup.R:69` (`if (by_mask == TRUE && methods::is(raster_mask, "SpatRaster"))`, otherwise the distance branch). Also `README.Rmd:242`, which passes `raster_mask = tenm_mask[1]`.

**What's wrong:** if `raster_mask` is missing or not a `SpatRaster`, the function quietly runs the distance-based method with `threshold`. `tenm_mask[1]` returns the value of cell 1 as a `data.frame`, so the README's headline time-specific thinning example (the step that keeps records from different years in the same pixel) actually ran distance thinning. The record counts quoted in the README come from the wrong method.

**Evidence:**

```r
class(tenm_mask[1])                               #> "data.frame"
nrow(clean_dup_by_date(abt, by_mask = TRUE, threshold = res(tenm_mask)[1],
                       raster_mask = tenm_mask[1])$temporal_df)   #> 40
nrow(clean_dup_by_date(abt, by_mask = TRUE, threshold = res(tenm_mask)[1],
                       raster_mask = tenm_mask)$temporal_df)      #> 50
```

**Fix:** `if (by_mask && !inherits(raster_mask, "SpatRaster")) stop("'raster_mask' must be a SpatRaster when by_mask = TRUE")`. Also fix `README.Rmd:242` and re-knit, because some README numbers will change.

---

### B9 (bug): distance thinning over-thins chains of nearby records

**Where:** `R/clean_dup.R:150-162`

**What's wrong:** the code builds the lower-triangular "within threshold" matrix and removes every record `j` that has any later neighbour `i`, even when `i` is itself removed. In a chain A–B–C (A–B and B–C closer than the threshold, A–C farther), it keeps only C, although {A, C} is a valid thinned set. On real clustered data this drops more records than intended.

**Evidence:**

```r
clean_dup(data.frame(id = 1:3, x = c(0, 0.6, 1.2), y = 0), "x", "y", threshold = 1)
#>   id   x y
#> 3  3 1.2 0      # points 1 and 3 are 1.2 apart, so greedy thinning keeps both
```

**Fix:** use greedy sequential thinning: walk the records and keep one only if it is not within the threshold of an already kept record. Use the sparse index from `sf::st_is_within_distance(x, dist = threshold)`, which also fixes P3. If the current behaviour is intended, document it.

---

### B10 (bug): points outside the mask are collapsed or dropped without warning

**Where:** `R/clean_dup.R:76-77` (`duplicated(cellids)` treats all `NA` cells as duplicates of each other) and `:141` (`dat2$cellid %in% keep` drops all `NA`).

**Evidence:** three points far outside the mask extent and one inside:

```r
out <- data.frame(id = 1:4, lon = c(-120, -121, -122, -97.3), lat = c(40, 41, 42, 18.7))
clean_dup(out, "lon", "lat", by_mask = TRUE, raster_mask = tenm_mask, n_ngbs = 0)  #> keeps ids 1, 4
clean_dup(out, "lon", "lat", by_mask = TRUE, raster_mask = tenm_mask, n_ngbs = 1)  #> keeps id 4 only
```

**Fix:** count `NA` cells, `warning()` with the count, and drop them consistently in both branches (or keep them with an option). This often signals a CRS mismatch (see R3).

---

### B11 (bug): `bg_by_date(buffer_distance = ...)` uses the `sample()` scalar trap

**Where:** `R/bg_by_date.R:100-102`. A similar call is `R/cells2samp.R:89`, which is reachable when the buffer yields a single cell.

**What's wrong:** when the buffer covers exactly one cell (`length(cell2sa) == 1`), `sample(cell2sa, 1)` samples from `1:cell2sa`, which returns an arbitrary cell anywhere in the raster. Overlapping buffers also return duplicate cells, because `cells()` is not `unique()`d.

**Evidence:** an occurrence at a cell centre with the default `buffer_distance = 1000` (m), on ~18 km cells:

```r
#> cells in 1 km buffer: 272
#> sample(cell2sa, 1) returned: 198  (value there: NA)
bg_by_date(one, buffer_ngbs = NULL, buffer_distance = 1000, n_bg = 10)
#> Error in tidyr::pivot_wider(): ! Must select at least one item.
```

**Fix:** `cell2sa <- unique(cell2sa); cell_ids <- cell2sa[sample.int(length(cell2sa), nsamples)]`. Also reconsider the default `buffer_distance = 1000`, which is metres for lon/lat rasters and smaller than one pixel of the example data.

---

### B12 (bug): `pROC()` ignores `rseed` and resets the global RNG

**Where:** `R/PartialROC.R:131` (`if (is.numeric(rseed)) set.seed(rseed)`) and `:193` (always called with `rseed = i`).

**What's wrong:** each bootstrap iteration calls `set.seed(i)` whatever the user's `rseed` argument. After `pROC()` returns, the global RNG is always in the state left by `set.seed(n_iter)`. Inside `ellipsoid_selection()` this means every model after the first starts `MASS::cov.mve` from the same seed. It also replaces any seed the user set before calling `pROC()`. CRAN policy discourages packages from changing the user's RNG state.

**Evidence:**

```r
set.seed(111); invisible(pROC(bgv, tv, n_iter = 20, rseed = FALSE)); a <- runif(1)
set.seed(999); invisible(pROC(bgv, tv, n_iter = 20, rseed = FALSE)); b <- runif(1)
identical(a, b)   #> TRUE   (0.2892624 both times)
```

**Fix:** pass `rseed = if (isTRUE(rseed)) i else NULL`, and when seeding, save and restore `.Random.seed` (for example with `withr::with_seed()`). Better still, draw all bootstrap indices once, as a `n_samp × n_iter` matrix, from the caller's stream.

---

### B13 (bug): `ellipsoid_omr()` finds objects in the global environment

**Where:** `R/ellipsoid_omr.R:182` (`if (exists("in_etest"))`).

**What's wrong:** `exists()` defaults to `inherits = TRUE`, and the namespace's parent chain reaches `.GlobalEnv`. If the user has any object called `in_etest`, a call without `env_test` takes the test branch and crashes.

**Evidence:**

```r
in_etest <- "anything"
ellipsoid_omr(tr, env_test = NULL, env_bg = bg, cf_level = 0.975, mve = FALSE)
#> Error: object 'occs_fail_test' not found
```

**Fix:** set a flag in the test block (`has_test <- TRUE`) and test that flag instead.

---

### B14 (bug): `tenm_selection()` documents parameters that do nothing

**Where:** `R/tenm_selection.R:109-112`

**What's wrong:** `RandomPercent`, `sub_sample` and `sub_sample_size` are documented as controlling the pROC bootstrap but never appear in the body. `pROC()` always runs with `boost_percent = 50`, `sub_sample = TRUE` and `sub_sample_size = 1000`. `E_percent` (fixed at 5) cannot be set at all. Together with B1, none of the pROC settings a user passes to `tenm_selection()` take effect.

**Evidence:** I searched the function body for each name: `RandomPercent: FALSE`, `sub_sample: FALSE`, `sub_sample_size: FALSE`.

**Fix:** pass them through `ellipsoid_selection()` → `ellipsoid_omr()` → `pROC(boost_percent = RandomPercent, sub_sample = ..., sub_sample_size = ..., E_percent = ...)`, or deprecate them.

---

### B15 (bug): several `predict()` arguments are broken

**Where:** `R/00_methods.R`

1. `:151, :177, :203`: `plot = TRUE` is hard-coded and `...` is also forwarded, so `predict(..., plot = FALSE)` errors with `formal argument "plot" matched by multiple actual arguments`. Users cannot turn off the `persp`/`rgl` plot, which opens a device for every projected layer.
2. `:133`: the error message uses `model_vars[idvars]` (the list of all models' variables) instead of `mod_vars[idvars]`, so it prints `Not valid variable names: NULL`.
3. `:193`: the default `layers_ext = NULL` reaches `gsub(NULL, ...)` and errors with `invalid 'pattern' argument`.

**Evidence:**

```r
predict(mod_sel, model_variables = c("bio_01","bio_12"), layers_path = dir70, layers_ext = ".tif$", plot = FALSE)
#> Caused by error in tenm::ellipsoid_projection(): ! formal argument "plot" matched by multiple actual arguments
predict(mod_sel, model_variables = c("bio_01","bio_99"), layers_path = dir70, layers_ext = ".tif$")
#> Error: Not valid variable names: NULL please provide valid variable names
predict(mod_sel, model_variables = c("bio_01","bio_12"), layers_path = dir70)
#> Caused by error in gsub(): ! invalid 'pattern' argument
```

**Fix:** add a `plot = FALSE` argument to `predict()` and forward it. Use `mod_vars[idvars]` in the message. Default `layers_ext` to the object's `layers_ext` (`object$layers_ext`).

---

### B16 (bug): `sp_temporal_data()` rejects tibbles and the documented `"mdy"` format

**Where:** `R/sp_temporal_data.R:60` and `:21-23,76`

**What's wrong:** `if (class(occs) %in% "data.frame")` has length 3 for a tibble, which is an error since R 4.2. `occs[, sp_date_var]` would also return a tibble rather than a vector. Separately, the documentation lists `"mdy"`, but `date_formats` contains `"myd"` and not `"mdy"`.

**Evidence:**

```r
sp_temporal_data(tibble::as_tibble(abronia), ...)
#> Error in if (class(occs) %in% classes_sp_data) { : the condition has length > 1
sp_temporal_data(occ, ..., occ_date_format = "mdy", ...)
#> occ_date_format should have one of these values y ym ymd my myd dmy
```

**Fix:** use `if (!inherits(occs, "data.frame")) stop(...)` and `occs <- as.data.frame(occs)`. Add `"mdy"` (`lubridate::mdy`), or correct the documentation.

---

### B17 (bug): exact duplicate records break `ex_by_date()`

**Where:** `R/ex_by_date.R:84-86`

**What's wrong:** `pivot_wider()` identifies rows by `(lon, lat, date, layer_dates, layers_path, cell_ids_year)`. Two identical records (common in GBIF downloads when `clean_dup_by_date()` is skipped) therefore produce list-columns with only a warning. `env_data`, `correlation_finder()` and `tenm_selection()` then fail later with unrelated errors.

**Evidence:**

```r
x <- abronia[c(1, 1, 2:6), ]; x$year <- 2014
ex <- ex_by_date(sp_temporal_data(x, ...), 0.7)
#> Warning: Values from `layer_val` are not uniquely identified; output will contain list-cols.
class(ex$temporal_df$bio_01)   #> "list"
```

**Fix:** add a row identifier (`tdf$.row <- seq_len(nrow(tdf))`) before extraction and use it as `id_cols`. Better, extract each year's stack in one `terra::extract()` call, which needs no reshaping (see P5).

---

### B18 (bug): `cov_center(vars = NULL)` (the default) fails with a misleading error

**Where:** `R/cov_center.R:53` (`data <- data[, vars]`)

**Evidence:**

```r
cov_center(etr, mve = TRUE, level = .975)
#> Error in MASS::cov.rob(...): 'x' is probably collinear
```

**Fix:** `if (!is.null(vars)) data <- data[, vars, drop = FALSE]`. Also validate that the data are numeric, have no `NA`, and have `nrow > ncol`, and use `tryCatch(solve(vari))` (`:80`) to give an informative "singular covariance" message.

---

### B19 (bug): `non_pred_test_ids` is `""` instead of `NA`

**Where:** `R/ellipsoid_omr.R:124`. The code tests `length(fails_train_ids) > 0`, which is always `TRUE` because that object is a length-1 string or `NA`, instead of `fails_test_ids`.

**Evidence:** with 0 of 5 test points omitted, `non_pred_test_ids` is `""`.

**Fix:** test `length(fails_test_ids) > 0`.

---

### B20 (bug): `clean_dup()` fails on common input shapes

**Where:** `R/clean_dup.R:65` (only `NA` longitudes are removed) and `:85-87` (the hard-coded `c(3, 1, 2, ...)` column order).

**Evidence:**

```r
x <- abronia[1:5, ]; x$decimalLatitude[2] <- NA
clean_dup(x, "decimalLongitude", "decimalLatitude", threshold = .1)
#> Error in st_as_sf.data.frame(...): missing values in coordinates not allowed
clean_dup(abronia[, 2:3], "decimalLongitude", "decimalLatitude", by_mask = TRUE, raster_mask = tenm_mask)
#> Error in `[.data.frame`(dat2, ids_nodup, c(3, 1, 2)): undefined columns selected
```

**Fix:** drop rows with `NA` in either coordinate (with a message). Subset the original `data` by `ids_nodup` instead of rebuilding it from `sf` and reordering by position. This also fixes the column reordering described in R6.

---

### B21 (bug): `proc_precision()` errors on degenerate test data

**Where:** `R/PartialROC.R:217` (`unique(sort(test_data))[2]` is `NA` when there are fewer than 2 unique values), and `:232` (`which(med == ".")` is empty when the number prints without a decimal point).

**Evidence:**

```r
pROC(bgv, rep(0.5, 10), n_iter = 5)
#> Error: missing value where TRUE/FALSE needed
```

Inside `ellipsoid_omr()` this is caught by `try()`, so the pROC columns silently become `NA` and the model drops to the bottom of the ranking.

**Fix:** compute precision numerically, for example `ndigits <- max(3, ceiling(-log10(min(diff(sort(unique(c(vals, test_data))))))))`, with a fallback when there is no spread. Or avoid rounding altogether by computing the ROC on sorted unique thresholds.

---

## B. Risks (silent or fragile behaviour)

### R1 (risk): occurrences without a layer for their date are silently re-dated

**Where:** `R/sp_temporal_data.R:69-75, 114-118`

**What's wrong:** each occurrence is matched to the layer folder with the nearest date, with no maximum gap and no message. Dates outside the range of available layers clamp to the first or last folder. Records whose dates fail to parse are dropped silently. A wrong `layers_date_format` or a non-existent directory gives the cryptic error `Result must be length 1, not 0`. Because `clean_dup_by_date()` groups by `layers_path`, not by date (`clean_dup_by_date.R:77`), two records from different years that both snap to the same layer are also deduplicated against each other. That runs against the package's premise of keeping same-pixel records from different years.

**Evidence:**

```r
occ <- data.frame(lon = -97.3, lat = 18.7, year = c(1945, 1939, 1890, 2050, 2016))
sp_temporal_data(occ, "lon", "lat", "year", "y", "y", d)$temporal_df[, c("year","layer_dates")]
#>   year layer_dates
#> 1 1945  1941-01-01
#> 2 1939  1939-01-01
#> 3 1890  1939-01-01
#> 4 2050  2016-01-01
#> 5 2016  2016-01-01
# Unparseable dates: input rows 5 -> output rows 3, no message.
# layers_date_format = "dmy":  "All formats failed to parse" warning, then "Result must be length 1, not 0".
# layers_by_date_dir = "/no/such/dir":  "Result must be length 1, not 0".
```

**Fix:** add `max_date_gap` (default 0, meaning exact match required) and `on_unmatched = c("drop", "nearest", "error")` arguments. Report how many records were matched exactly, approximately, or dropped. Check `dir.exists()` and `all(is.na(layers_dates))` up front with explicit messages. Keep the original date in the output (it is already kept) and add a `date_gap` column.

---

### R2 (risk): `predict()` does not project the model that was selected

**Where:** `R/00_methods.R:113-143`

**What's wrong:** `tenm_selection()` stores only the metrics table, not the fitted centroid and covariance. `predict()` refits the model with its own `mve = TRUE, level = 0.975` defaults. If selection used `ellipsoid_level = 0.95` or `mve = FALSE`, the projection silently uses a different ellipsoid. Because `MASS::cov.mve` is stochastic, even matching settings give a slightly different ellipsoid from the one that was evaluated.

**Fix:** store `ellipsoid_level`, `mve` and the seed in the `sp.temporal.selection` object, and use them as `predict()` defaults. Ideally, also store `centroid` and `covariance` for the top models.

---

### R3 (risk): CRS, extent, resolution and layer completeness are never validated

**Where:** `sp_temporal_data`, `ex_by_date`, `bg_by_date`, `clean_dup`, `predict`.

**What's wrong:**

- The occurrence CRS is assumed to equal the raster CRS. `clean_dup` just stamps the mask's CRS on the coordinates (`clean_dup.R:70-71`).
- In distance mode, `threshold` is in planar coordinate units, which means degrees for lon/lat. A degree is a different distance at different latitudes (`clean_dup.R:148-149`).
- With projected layers and lon/lat occurrences, every `cellFromXY` returns `NA`. That leads to B10 (records silently dropped) and to all-`NA` environmental values.
- If a year folder lacks one variable, `pivot_wider` fills that column with `NA` for that year without a message.
- Grids that differ across years give wrong values (B6).

**Fix:** add one `validate_layers()` helper, run in `sp_temporal_data()`, that checks for every folder:
- the same set of layer names;
- `terra::compareGeom()` against the first folder;
- that the CRS is defined.

Accept an `occ_crs` argument and `terra::project()` the points if needed. In distance mode, use `sf::st_is_within_distance()` on lon/lat points with s2, so the threshold is in metres.

---

### R4 (risk): ranking hides failures and missing metrics

**Where:** `R/ellipsoid_selection.R:251-254, 284`, `R/ellipsoid_omr.R:69-77`

**What's wrong:**
- A missing `om_rate_test` (for example, no test records for those variables) is replaced by 0, so the model looks perfect on test data.
- Models whose fit fails (singular covariance, too few points for `cov.mve`) return `NULL` and disappear from the table with no record.
- When no model meets `omr_criteria`, the function returns early without the `rank_omr_aucratio` column, so the output schema depends on the data.

**Fix:** keep `NA` and rank `NA` last. Return a row with `status = "fit_failed"` and the error message for failed models. Always return the same columns.

---

### R5 (risk): background prevalence includes the training presences

**Where:** `R/ellipsoid_omr.R:163` (`table(c(in_ebg$in_Ellipsoid, in_e$in_Ellipsoid))`)

**What's wrong:** prevalence is estimated over background plus training presences, and about 97.5% of the presences are inside by construction. This inflates `bg_prevalence`, which is also the null probability in the binomial test. If this is intentional (similar to Maxent's "add samples to background"), please document it. Otherwise use only `in_ebg`.

---

### R6 (risk): columns addressed by position

**Where:** `ex_by_date.R:77` (`time_obs[, c(1:6)]`) and `:111` (`years_env[, -(1:6)]`); `tdf2swd.R:41` (`-c(4, 5, 6, ncol)`) and `:48` (`df_base[, c(2, 3)]`); `bg_by_date.R:147-148` (`bg_env[-(1:2)]`, `bg_env[[3]]`).

**What's wrong:** these assume a fixed column layout. `clean_dup(by_mask = TRUE)` already moves the date column to the front (`clean_dup.R:85`). As a result, `tdf2swd()` emits `sp_name, year, lon, lat, ...` instead of the SWD layout `sp_name, lon, lat`. `env_data` also carries the `trian_test` label as its last column, which is why the README and examples must call `correlation_finder(abex$env_data[, -ncol(abex$env_data)])`. In addition, `tdf2swd()` derives `year` from the folder name (`tdf2swd.R:44`), which gives `NA` with a warning for any non-`"y"` folder name (182 `NA`s in my test with a `1939-06` folder).

**Fix:** select by name (`this_species$lon_lat_vars`, `sp_date_var`, a stored `env_vars` vector). Exclude `trian_test` from `env_data`. Use `layer_dates`, not the folder name, for `year`.

---

### R7 (risk): code depends on terra internals through `eval(parse())`

**Where:** `ex_by_date.R:73-74` and `bg_by_date.R:130-131` (`env_layers@<slot>$get_sourcenames()`), `PartialROC.R:63-66` (`@.xData$range_min`/`range_max`), `metaras.R:13-15` (`@.xData$names`).

**What's wrong:** these reach into terra's C++ module fields. The S4 slot that holds the module is an undocumented internal, and it has already changed name: it is `cpp` in terra 1.7-65 and `pntr` in terra 1.9-50 (both checked). The code still works on both versions only because it takes `names(attributes(x))[1]`. Any change to the module's method or field names will break the package. In `pROC()`, `range_min` is also `NA` when min/max have not been computed, which makes `if (ck_min == ck_max)` error.

**Fix:** use the public API: `names(env_layers)` (what `metaras()` returns, which I checked) or `terra::sources()`, and `terra::minmax(x, compute = TRUE)` or `terra::global(x, "range")`.

---

### R8 (risk): the parallel branch of `ellipsoid_selection()` changes global state and is not reproducible

**Where:** `R/ellipsoid_selection.R:132, 146-151, 175-177, 223`

**What's wrong:**

- `options(future.rng.onMisuse = "ignore")` and `options(future.globals.maxSize = 8.9e9)` are set and never restored. CRAN policy forbids changing options without resetting them.
- `plan(tweak(multisession, workers = n_cores))` overrides the user's plan. `future` guidance is that packages should not set plans; they should respect the user's.
- `availableCores() - 1` can be 0 on a single-core container.
- Futures are created without `seed = TRUE`, and RNG-misuse warnings are suppressed, so `cov.mve` and pROC use unsound parallel RNG streams.

**Evidence:**

```r
before: future.rng.onMisuse = NULL, future.globals.maxSize = NULL
after : future.rng.onMisuse = "ignore", future.globals.maxSize = 8.9e+09
set.seed(10); ra <- ellipsoid_selection(..., parallel = TRUE, ncores = 2)
set.seed(10); rb <- ellipsoid_selection(..., parallel = TRUE, ncores = 2)
identical(ra, rb)   #> FALSE   (serial: TRUE)
```

**Fix:** use `furrr::future_map(..., .options = furrr_options(seed = TRUE))` (or `future.apply` with `future.seed = TRUE`) under whatever plan the user has set. Use `withr::local_options()` if an option is truly needed. Use `max(1, availableCores() - 1)`.

---

## C. Performance

### P1 (improvement): `bg_by_date()` spends most of its time unpacking list-columns row by row

**Where:** `R/bg_by_date.R:139-159`

**What's wrong:** `pivot_wider(values_fn = list)` creates list-columns. The code then rebuilds the table with `future_map_dfr` over every row, and for each row a `map_dfc` over every variable. That is about 9 ms per background point. The default `n_bg = 50000` therefore takes roughly 7-8 minutes for this step alone.

**Evidence:** `Rprof` with `n_bg = 4000` shows `purrr::map_dfc` at 85% of total time. Timings: `n_bg = 2000` took 22.9 s and `n_bg = 8000` took 73.2 s. On a synthetic table of the same shape (50,000 rows × 19 list-columns), `tidyr::unnest(bg_env, all_of(vars))` takes 0.74 s.

**Fix:** extract per year with one call, `terra::extract(terra::rast(files_of_year), cells)` (returns a wide `data.frame`), then `rbind` the years. This removes `pivot_wider`, the list-columns and the per-row loop. Remember to take `unique()` of the sampled cells first.

### P2 (improvement): pROC dominates model selection

**Where:** `ellipsoid_omr.R:193`, `PartialROC.R:186-194`

**What's wrong:** `pROC()` runs for every candidate model, and each iteration builds a `data.frame` that is `map_df`-bound. It takes 83% of selection time: 0.21 s per model with `proc = TRUE` against 0.04 s without, and that is at the accidental 100 iterations (B1). At the documented 1000 iterations it is about 2 s per model. With 19 variables and `nvars_to_fit = 2:7` (94,164 models) that is days of CPU.

**Fix:**
- Use two-stage selection: compute omission rates for all models (cheap), then run pROC only for models that pass `omr_criteria`.
- Vectorize `pROC`: draw an `n_samp × n_iter` index matrix once, compute sensitivities with `findInterval()` on sorted thresholds, and return numeric vectors instead of `data.frame`s.
- Cache `qchisq` and the background `mahalanobis` computations per variable combination.

### P3 (improvement): `clean_dup()` distance mode needs quadratic memory

**Where:** `R/clean_dup.R:148-150` (`st_buffer`, then `st_intersects(sparse = FALSE)`)

**Evidence:** n = 2,000 used 225 MB in 2.3 s; n = 4,000 used 393 MB in 3.9 s; n = 8,000 used 1,245 MB in 14.1 s. Extrapolating, n = 30,000 would need about 17 GB.

**Fix:** use `sf::st_is_within_distance(x, dist = threshold)` (sparse, spatially indexed; no buffers needed), combined with the greedy fix from B9.

### P4 (improvement): `ellipsoid_projection()` loads entire rasters into memory

**Where:** `R/ellipsoid_projection.R:64-96`

**What's wrong:** `envlayers[[x]][]` reads each layer fully and `map_dfc` builds a `data.frame` column by column. `mahalanobis()` then makes further copies. For 1 km continental or global layers this will not fit in RAM.

**Fix:** use `terra::app(envlayers, fun = function(v) exp(-0.5 * mahalanobis(v, centroid, covar)))` or `terra::predict()` with a small model object. Both process in blocks and can write to disk (`filename =`).

### P5 (improvement): redundant raster I/O

**Where:** `ex_by_date.R:56-81` and `bg_by_date.R:93, 118-137`

**What's wrong:**
- `ex_by_date` opens every file separately, filters the whole `tdf` once per file (years × variables tasks), and reshapes the result. This took 9 s for 50 records × 32 years × 19 variables.
- `bg_by_date`'s `buffer_distance` branch calls `terra::rast(ddL[[z]]$layer_path)` with one identical path per occurrence, building an n-layer raster.

**Fix:** use one task per year: `terra::extract(terra::rast(files_of_year), xy)`. Use `terra::rast(ddL[[z]]$layer_path[1])`.

---

## D. Dependencies and CRAN-check risks

- **D1 (good):** no deprecated spatial dependencies. `raster`, `sp`, `rgdal`, `rgeos` and `maptools` are not used; the only reference is the commented-out line `clean_dup.R:68`. All spatial work uses `terra` and `sf`.
- **D2 (improvement):** global options are changed and not restored (R8). CRAN policy says "functions ... should not change the user's options ... without restoring them". `pROC()` also re-seeds the RNG (B12).
- **D3 (improvement):** `furrr`'s `.progress` argument is documented as "will be deprecated and removed in a future version of furrr in favor of using the progressr package". It is hard-coded to `TRUE` in `clean_dup_by_date.R:90`, `ex_by_date.R:81` and `bg_by_date.R:137,157`, ignoring the function's own `progress` argument, and floods logs and knitted output with progress bars. Switch to `progressr` and respect `progress`.
- **D4 (improvement):** `rgl` is in `Imports` but is only used for optional 3D plots (`ellipsoid_projection`, `plot_ellipsoid`). It pulls in X11/OpenGL system libraries and is a common source of install failures on headless servers and HPC systems. Move it to `Suggests` and guard with `requireNamespace("rgl", quietly = TRUE)`.
- **D5 (improvement):** unversioned dependencies on features that need recent versions:
  - `purrr::list_rbind()` and `purrr::map(.progress =)` need purrr ≥ 1.0.0;
  - `dplyr::summarise(.by =)`, if adopted, needs dplyr ≥ 1.1.0.

  Declare `purrr (>= 1.0.0)`. `MASS` could also move to a lighter robust-covariance option, but that is optional.
- **D6 (improvement):** repository and test hygiene:
  - `tests/testthat/Rplots.pdf` is committed and ships in the built tarball; I checked with `tar tzf`.
  - `.RData`, `.Rhistory` and `.DS_Store` are committed. `R CMD build` excludes them, but they should be removed and added to `.gitignore`.
  - The `tenm_selection` example uses `n_cores = 20`. CRAN allows at most 2 cores in examples; this is only saved by the `availableCores() - 1` cap, and only when `_R_CHECK_LIMIT_CORES_` is honoured.

**R CMD check results** (`R CMD check --as-cran --no-manual --ignore-vignettes`, R 4.3.3, `_R_CHECK_CRAN_INCOMING_=false`):

- **0 errors, and all tests pass** (`testthat.R` took 70 s CPU / 137 s elapsed).
- **Examples:** 12 s; with `--run-donttest`, 3 min CPU but **11 min elapsed**, mostly from `bg_by_date(n_bg = 50000)` (P1), which several examples call. CRAN's `--run-donttest` checks will be slow. Use `n_bg = 1000` in examples.
- **NOTE (real):** installed size is 5.6 MB, of which `extdata` is 5.0 MB. Consider fewer or coarser example layers, or compressing the GeoTIFFs (`gdal = "COMPRESS=DEFLATE"`).
- **WARNING and NOTE (environment only, not package problems):** "cannot set locale en_US.UTF-8" and "unable to verify current time" both come from the sandbox container.


---

## E. Maintainability

- **M1: `ellipsoid_selection()` (~200 lines) duplicates its core loop.** The serial (`:226-246`) and parallel (`:189-212`) branches each reimplement "subset, `na.omit`, `ellipsoid_omr`". That duplication is how B1 got into both copies. Extract `fit_one(vars)`, build the list of combinations once, and call a single `furrr::future_map()` (which runs sequentially under the default plan). The NA-padding of `combn` matrices and the manual chunking via `kkk`/`pasos` then go away.
- **M2: copy-pasted date parsing.** `sp_temporal_data.R:80-112` has two identical six-way `if`/`else` chains. Use one lookup table, `parsers <- list(y = function(x) lubridate::ymd(x, truncated = 2L), ym = lubridate::ym, ...)`, and `parsers[[fmt]](x)`.
- **M3: `eval(parse())` used for ordinary variables.** `cov_center.R:95-101` uses `assign()` + `eval(parse())` to build a 2-row matrix: `rbind(vec_1 = centroid + v * s, vec_2 = centroid - v * s)`. R7 lists the uses that reach into terra internals.
- **M4: repeated counting logic.** `ellipsoid_omr.R` uses the `table()`/`which(names(...) %in% "1")` pattern three times to count 0/1 values. Use `sum(x == 1, na.rm = TRUE)` / `sum(x == 0, na.rm = TRUE)`, and decide explicitly how `NA` rows enter the denominator (today they count in `nrow()` but not as omissions, which lowers the omission rate).
- **M5: one source of truth for the layer inventory.** `ex_by_date`, `bg_by_date` and `predict` each rediscover files with `list.files()` and regex juggling (`ex_by_date.R:56-64` strips file names with an unescaped regex built from basenames). Build a validated `layers_index` (year × variable × path) once in `sp_temporal_data()`, store it in the object, and use it everywhere. This also fixes B4, B5, B6 and R3 structurally.
- **M6: API consistency.**
  - The same concepts have different names in different functions:
    - `n_cores` vs `ncores`;
    - `NoOfIteration`/`RandomPercent` vs `proc_iter`/`boost_percent`;
    - `ellipsoid_level` vs `level` vs `cf_level`.
  - The column name `trian_test` (a typo) is part of the public output.
  - `env_data` includes the label column (R6).
  - Classes are inconsistent:
    - `bg_by_date()` returns `class = c("sp.temporal.env", "sp.temporal.bg")`, but `env_data` is dropped;
    - `tenm_selection()` returns only `"sp.temporal.selection"`, although the docs say it inherits `sp.temporal.bg` (`inherits(x, "sp.temporal.bg")` is `FALSE`);
    - the S4 classes in `00_classes.R` (with `validity = TRUE` stubs) are never instantiated, because every object is an S3 list.

  Pick S3, give each class a constructor and validator, and use `stopifnot(inherits(...))` consistently (`ex_by_date()` has no class check).
- **M7: tests do not catch regressions.** `tests/testthat/test-tenm.R` checks only classes (`expect_s3_class`, `expect_match(class(...))`), so none of the bugs above would fail a test. Suggested regression tests, each a few lines with the bundled data:
  - B1: `mve = FALSE` is honoured.
  - B2: `predict()` is invariant to variable order.
  - B3: `pval_bin` equals `binom.test(..., "greater")`.
  - B4/B5: output is invariant to row order and folder naming.
  - B8: an invalid `raster_mask` errors.
  - B9: chain thinning keeps 2 of 3.
  - B12: RNG state is unchanged after `pROC(rseed = FALSE)`.

---

## Suggested order of work

1. **One-line fixes with large impact:** B1, B3, B13, B19, B15.2 and B8 (including re-knitting the README).
2. **Projection correctness:** B2, then R2 (store model parameters in the selection object).
3. **Alignment bugs:** B4, B5, B6 and R3, via a validated layer index (M5).
4. **Curation:** B7, B9, B10 and R1 (explicit handling of unmatched dates).
5. **Performance:** P1, then P2 (two-stage pROC), which makes realistic model sets feasible.
6. **CRAN hygiene and API cleanup:** section D and M6, with the regression tests from M7.
