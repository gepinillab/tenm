# tenm 0.5.1.9002

## Bug fixes

* `ellipsoid_omr()` no longer takes the test-data branch when an object
  called `in_etest` exists in the global environment. Before, a call
  without `env_test` failed with "object 'occs_fail_test' not found" (#14).

* `tenm_selection()` now passes `RandomPercent`, `sub_sample` and
  `sub_sample_size` to `pROC()`, and gains an `E_percent` argument. Before,
  they were ignored and `pROC()` always ran with `boost_percent = 50`,
  `sub_sample = TRUE`, `sub_sample_size = 1000` and `E_percent = 5`.
  `ellipsoid_selection()` and `ellipsoid_omr()` gain `E_percent`,
  `boost_percent`, `sub_sample` and `sub_sample_size` arguments with those
  values as defaults, so results with default arguments are unchanged (#15).

* `predict()` gains a `plot` argument (default `FALSE`). Before, `plot = TRUE`
  was hard-coded, so `predict(..., plot = FALSE)` failed and every projected
  layer opened a plot. The invalid-variable error now names the invalid
  variables instead of printing `NULL`, and `layers_ext = NULL` now falls
  back to the object's `layers_ext` instead of failing in `gsub()` (#16).

* `ellipsoid_omr()`: `non_pred_test_ids` is now `NA` when no test point is
  omitted. Before, it was `""` (#20).

* `pROC()` no longer fails when all test values are equal, or when the
  values have no decimal part (e.g. 0-100 suitabilities). Inside
  `ellipsoid_omr()` these errors made the partial ROC columns silently `NA`.
  Results for inputs that worked before are unchanged (#22).

* `tenm_selection()` now stores `mve` and `ellipsoid_level` in the
  `sp.temporal.selection` object, and `predict()` uses them when `mve` and
  `level` are not supplied. Before, `predict()` always refit with
  `mve = TRUE, level = 0.975`, whatever settings were used in selection (#24).

* `ellipsoid_selection()` keeps a missing test omission rate (no test
  records for a variable combination) as `NA` instead of 0, and ranks those
  models after the models with test data. Models whose ellipsoid cannot be
  fitted (e.g. singular covariance) are now listed last with
  `status = "fit_failed"` and the error message in `fit_error`, instead of
  being dropped; `ellipsoid_omr()` returns such a row instead of `NULL`.
  `rank_omr_aucratio` is now present (with `proc = TRUE`) even when no model
  meets `omr_criteria`, and every table has `status` and `fit_error`
  columns (#26).

* `ellipsoid_omr()`: `bg_prevalence`, which is also the null probability of
  the binomial test (`pval_bin`), is now the proportion of background points
  inside the ellipsoid. Before, the training presences were counted as
  background, which inflated `bg_prevalence` and `pval_bin` (#27).

# tenm 0.5.1.9001

## Bug fixes

* `ellipsoid_selection()` (and so `tenm_selection()`) now passes `mve` and
  `proc_iter` to `ellipsoid_omr()` by name. Before, `proc_iter` bound to
  `mve`, so `mve = FALSE` was ignored and the partial ROC always ran 100
  iterations whatever `proc_iter`/`NoOfIteration` was (#2).

* `predict()` and `ellipsoid_projection()` now match layers to model
  variables by name. Before, layers were used in file or stack order, so
  any non-alphabetical variable order gave a wrong map without an error.
  A variable with no layer is now an informative error (#3).

* `ellipsoid_omr()`: `pval_bin` is now the one-sided binomial p-value
  P(X >= k) (as in `binom.test(..., alternative = "greater")`). Before it
  was P(X > k), which understated the p-value and was 0 whenever all test
  points were inside the ellipsoid (#4).

* `bg_by_date()` now aligns template layers, per-date sample sizes and
  sampled cells by date folder instead of by row position. Before, unsorted
  `temporal_df` rows or folder names that do not sort in date order
  (e.g. `"my"` folders) gave background from the wrong date and swapped
  sample sizes (#5).

* `ex_by_date()` now assigns train/test labels within each date. Before,
  labels landed on other dates' rows when folder order differed from date
  order (#6).

* `ex_by_date()` now stops with an informative error if the environmental
  layers do not all share one grid (extent, resolution and CRS). Before,
  the first date's cell IDs were reused for every date, so background
  values came from the wrong pixels (#7).

* `bg_by_date(buffer_distance = ...)` and `cells2samp()` no longer return an
  arbitrary cell when the buffer holds a single cell, and `bg_by_date()` no
  longer fails on overlapping buffers (#12).

* `pROC()` now honours `rseed`: with `rseed = FALSE` the bootstrap draws from
  the caller's random number stream, and with `rseed = TRUE` iterations are
  seeded for reproducibility and the caller's RNG state is restored on exit.
  Before, every call reset the global RNG to `set.seed(n_iter)` (#13).

# tenm 0.5.1

* Improvements to documentation

# tenm 0.5.0

* Initial CRAN submission.

* Time-specific niche modeling (TENM) is a novel approach that allows 
  calibrating niche models with high temporal resolution spatial information, 
  which aims to reduce niche estimation biases. 
  What makes the tenm package stand out is its distinctive capability to 
  calibrate models by incorporating specific temporal information, 
  whether on a yearly, monthly, or even daily basis. This feature distinguishes
  it from traditional models that rely on averaged temporal data. 
  
  Some of the package functions are:
  - Time-specific spatial data thinning: data cleaning considering the 
    temporal dimension of records. 
  - Time-specific environmental data extraction: extract environmental data 
    from variables considering the temporal dimension of data.
  - Time-specific background generation: generate background points for 
    the modeling process by considering the temporal dimension of the 
    occurrence and environmental data: The number of background points for
    each year is proportional to the number of occurrences for each year
    of observation. 
  - Exporting time-specific information as Samples With Data format: 
    export the time-specific data to the Samples With Data format table. 
    This function allows users to use other modeling algorithms such 
    as MaxEnt and GLMs.
  - Time-specific model calibration: calibrate time-specific niche models 
    using minimum volume ellipsoids. It fits numerous models based on a 
    combination of user-set parameters, including different combinations of
    environmental variables.
  - Time-specific model selection: select n models from all the fitted 
    models using statistical and model performances as model 
    selection criteria.
  - Projecting time-specific niche models: project one or more of the 
    selected models onto both the environmental and/or geographical space. 
