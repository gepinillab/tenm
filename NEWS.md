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
