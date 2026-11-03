# Regression tests for confirmed bugs (see gepinillab/tenm issues).

testthat::test_that("#2: ellipsoid_selection() honours mve and proc_iter", {
  set.seed(42)
  env_train <- data.frame(bio_01 = c(rnorm(50, 20, 2), rnorm(10, 30, 6)),
                          bio_12 = c(rnorm(50, 1000, 100), rnorm(10, 1500, 300)))
  env_test <- data.frame(bio_01 = rnorm(20, 21, 3),
                         bio_12 = rnorm(20, 1050, 150))
  env_bg <- data.frame(bio_01 = runif(3000, 5, 45),
                       bio_12 = runif(3000, 200, 2500))
  set.seed(1)
  direct <- tenm::ellipsoid_omr(env_data = env_train, env_test = env_test,
                                env_bg = env_bg, cf_level = 0.975,
                                mve = FALSE, proc = TRUE, proc_iter = 7,
                                rseed = TRUE)
  set.seed(1)
  utils::capture.output(
    sel <- tenm::ellipsoid_selection(env_train = env_train,
                                     env_test = env_test,
                                     env_vars = c("bio_01", "bio_12"),
                                     nvarstest = 2, level = 0.975,
                                     mve = FALSE, env_bg = env_bg,
                                     omr_criteria = 0.1, parallel = FALSE,
                                     proc = TRUE, proc_iter = 7,
                                     rseed = TRUE))
  cols <- c("om_rate_train", "om_rate_test", "bg_prevalence", "pval_bin",
            "pval_proc", "env_bg_paucratio", "env_bg_auc")
  testthat::expect_equal(unlist(sel[1, cols]), unlist(direct[1, cols]))
})

testthat::test_that("#3: predict() matches layers to model variables by name", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  dir70 <- system.file("extdata/bio_1970_2000", package = "tenm")
  bio_01 <- terra::rast(file.path(dir70, "bio_01.tif"))
  bio_12 <- terra::rast(file.path(dir70, "bio_12.tif"))
  names(bio_01) <- "bio_01"
  names(bio_12) <- "bio_12"
  set.seed(3)
  cells <- which(!is.na(terra::values(bio_01)[, 1]) &
                   !is.na(terra::values(bio_12)[, 1]))
  cells <- sample(cells, 40)
  train <- data.frame(bio_01 = terra::values(bio_01)[cells, 1],
                      bio_12 = terra::values(bio_12)[cells, 1],
                      trian_test = "Train")
  mod_sel <- list(temporal_df = train,
                  env_bg = train[, c("bio_01", "bio_12")],
                  mods_table = data.frame(fitted_vars = "bio_01,bio_12"))
  class(mod_sel) <- "sp.temporal.selection"
  utils::capture.output({
    p1 <- predict(mod_sel, model_variables = c("bio_01", "bio_12"),
                  layers_path = dir70, layers_ext = ".tif$", mve = FALSE)
    p2 <- predict(mod_sel, model_variables = c("bio_12", "bio_01"),
                  layers_path = dir70, layers_ext = ".tif$", mve = FALSE)
    p3 <- predict(mod_sel, model_variables = c("bio_01", "bio_12"),
                  layers = c(bio_12, bio_01), mve = FALSE)
    p4 <- predict(mod_sel, model_variables = c("bio_01", "bio_12"),
                  layers = list(c(bio_12, bio_01)), mve = FALSE)[[1]]
  })
  mod <- tenm::cov_center(train[, c("bio_01", "bio_12")], mve = FALSE,
                          level = 0.975, vars = c("bio_01", "bio_12"))
  vals <- cbind(terra::values(bio_01)[, 1], terra::values(bio_12)[, 1])
  ok <- stats::complete.cases(vals)
  ref <- rep(NA_real_, nrow(vals))
  ref[ok] <- exp(-0.5 * stats::mahalanobis(vals[ok, ], mod$centroid,
                                           mod$covariance))
  for (p in list(p1, p2, p3, p4)) {
    testthat::expect_equal(terra::values(p)[, 1], ref)
  }
  # A variable missing from the layers is an informative error
  testthat::expect_error(
    tenm::ellipsoid_projection(envlayers = c(bio_01, bio_01),
                               centroid = c(bio_01 = 1, bio_05 = 1, bio_12 = 1),
                               covar = diag(3), plot = FALSE),
    "bio_05")
})

testthat::test_that("#4: pval_bin is the one-sided binomial P(X >= k)", {
  set.seed(4)
  env_train <- data.frame(bio_01 = rnorm(60, 20, 3),
                          bio_12 = rnorm(60, 1000, 150))
  env_bg <- data.frame(bio_01 = runif(2000, 5, 40),
                       bio_12 = runif(2000, 300, 2000))
  tests <- list(
    all_inside = data.frame(bio_01 = c(19, 20, 21, 20, 20.5),
                            bio_12 = c(1000, 980, 1020, 1050, 990)),
    some_outside = data.frame(bio_01 = c(rnorm(9, 20, 2), 35, 6, 38),
                              bio_12 = c(rnorm(9, 1000, 100), 1900, 350, 400)))
  for (env_test in tests) {
    r <- tenm::ellipsoid_omr(env_data = env_train, env_test = env_test,
                             env_bg = env_bg, cf_level = 0.975, mve = FALSE)
    n <- nrow(env_test)
    k <- round(n * (1 - r$om_rate_test))
    ref <- stats::binom.test(k, n, p = r$bg_prevalence,
                             alternative = "greater")$p.value
    testthat::expect_equal(r$pval_bin, ref)
  }
})

testthat::test_that("#5: bg_by_date() is keyed by date folder, not row position", {
  future::plan("sequential")
  occ_xy <- data.frame(x = c(5.5, 24.5, 20.5, 24.5),
                       y = c(5.5, 24.5, 24.5, 20.5),
                       year_index = c(1, 2, 2, 2))
  scenarios <- list(list(folders = c("12-1939", "01-1940"), fmt = "my",
                         reverse = FALSE),
                    list(folders = c("1939", "1940"), fmt = "y",
                         reverse = TRUE))
  for (sc in scenarios) {
    lay <- make_year_layers(sc$folders)
    occs <- data.frame(lon = occ_xy$x, lat = occ_xy$y,
                       date = sc$folders[occ_xy$year_index])
    set.seed(5)
    st <- tenm::sp_temporal_data(occs = occs, longitude = "lon",
                                 latitude = "lat", sp_date_var = "date",
                                 occ_date_format = sc$fmt,
                                 layers_date_format = sc$fmt,
                                 layers_by_date_dir = lay$root,
                                 layers_ext = "*.tif$")
    ex <- tenm::ex_by_date(st, train_prop = 0.7)
    if (sc$reverse) {
      ex$temporal_df <- ex$temporal_df[rev(seq_len(nrow(ex$temporal_df))), ]
    }
    bg <- tenm::bg_by_date(ex, buffer_ngbs = 2, n_bg = 40)$env_bg
    year_index <- match(basename(bg$ID_YEAR), sc$folders)
    # sample size proportional to each year's occurrences (1 and 3 of 4)
    testthat::expect_equal(as.vector(table(factor(year_index, 1:2))),
                           c(10, 30))
    # each background cell lies within 2 cells of its own year's occurrences
    cell <- terra::cellFromXY(lay$template, as.matrix(bg[, c("lon", "lat")]))
    rc <- terra::rowColFromCell(lay$template, cell)
    occ_rc <- terra::rowColFromCell(
      lay$template,
      terra::cellFromXY(lay$template, as.matrix(occ_xy[, c("x", "y")])))
    in_buffer <- vapply(seq_len(nrow(bg)), function(j) {
      own <- occ_rc[occ_xy$year_index == year_index[j], , drop = FALSE]
      any(pmax(abs(own[, 1] - rc[j, 1]), abs(own[, 2] - rc[j, 2])) <= 2)
    }, logical(1))
    testthat::expect_true(all(in_buffer))
    # values come from the same year's layers at the reported coordinates
    testthat::expect_equal(bg$bio_01, year_index * 10000 + cell)
  }
})

testthat::test_that("#6: ex_by_date() assigns train/test labels within each date", {
  future::plan("sequential")
  folders <- c("12-1939", "01-1940", "03-1939")
  n_occ <- c(1, 2, 4)
  lay <- make_year_layers(folders)
  occs <- data.frame(lon = c(5.5, 10.5, 12.5, 20.5, 22.5, 24.5, 26.5),
                     lat = c(5.5, 10.5, 12.5, 20.5, 22.5, 24.5, 26.5),
                     date = rep(folders, n_occ))
  st <- tenm::sp_temporal_data(occs = occs, longitude = "lon",
                               latitude = "lat", sp_date_var = "date",
                               occ_date_format = "my",
                               layers_date_format = "my",
                               layers_by_date_dir = lay$root,
                               layers_ext = "*.tif$")
  for (seed in 1:3) {
    set.seed(seed)
    tdf <- tenm::ex_by_date(st, train_prop = 0.7)$temporal_df
    folder <- basename(tdf$layers_path)
    n_train <- vapply(folders, function(f) {
      sum(tdf$trian_test[folder == f] == "Train")
    }, numeric(1))
    testthat::expect_equal(unname(n_train), c(1, 1, 3))
    testthat::expect_equal(as.vector(table(factor(folder, folders))), n_occ)
  }
})

testthat::test_that("#7: ex_by_date() rejects layers that do not share one grid", {
  future::plan("sequential")
  occs <- data.frame(lon = c(10.5, 15.5, 20.5, 12.5),
                     lat = c(10.5, 15.5, 20.5, 22.5),
                     date = c("1939", "1939", "1940", "1940"))
  make_st <- function(lay) {
    tenm::sp_temporal_data(occs = occs, longitude = "lon", latitude = "lat",
                           sp_date_var = "date", occ_date_format = "y",
                           layers_date_format = "y",
                           layers_by_date_dir = lay$root,
                           layers_ext = "*.tif$")
  }
  # 1940 layers cropped by 3 columns: same CRS and resolution, other extent
  cropped <- make_year_layers(c("1939", "1940"), crop_folder = "1940")
  testthat::expect_error(tenm::ex_by_date(make_st(cropped), train_prop = 0.7),
                         "same extent, resolution and CRS")
  # Matching grids still work
  same <- make_year_layers(c("1939", "1940"))
  testthat::expect_s3_class(tenm::ex_by_date(make_st(same), train_prop = 0.7),
                            "sp.temporal.env")
})

testthat::test_that("#12: buffer_distance and cells2samp() handle single and repeated cells", {
  future::plan("sequential")
  lay <- make_year_layers("1939")
  bg_for <- function(occs, buffer_distance, n_bg) {
    st <- tenm::sp_temporal_data(occs = occs, longitude = "lon",
                                 latitude = "lat", sp_date_var = "date",
                                 occ_date_format = "y",
                                 layers_date_format = "y",
                                 layers_by_date_dir = lay$root,
                                 layers_ext = "*.tif$")
    ex <- tenm::ex_by_date(st, train_prop = 0.7)
    tenm::bg_by_date(ex, buffer_ngbs = NULL,
                     buffer_distance = buffer_distance, n_bg = n_bg)$env_bg
  }
  # A 1 km buffer around a cell centre covers only that cell
  occ <- data.frame(lon = 15.5, lat = 15.5, date = "1939")
  cell <- terra::cellFromXY(lay$template, cbind(15.5, 15.5))
  for (seed in 1:3) {
    set.seed(seed)
    bg <- bg_for(occ, buffer_distance = 1000, n_bg = 10)
    testthat::expect_equal(
      terra::cellFromXY(lay$template, as.matrix(bg[, c("lon", "lat")])), cell)
    testthat::expect_equal(bg$bio_01, 10000 + cell)
  }
  # Overlapping buffers do not yield duplicated background cells
  occ <- data.frame(lon = c(10.5, 11.5, 12.5), lat = 10.5, date = "1939")
  buf <- terra::buffer(terra::vect(as.matrix(occ[, 1:2]), crs = "EPSG:4326"),
                       width = 150000)
  n_cells <- length(unique(terra::cells(lay$template, buf)[, 2]))
  set.seed(1)
  bg <- bg_for(occ, buffer_distance = 150000, n_bg = 1000)
  testthat::expect_equal(nrow(bg), n_cells)
  testthat::expect_false(any(duplicated(bg[, c("lon", "lat")])))
  # cells2samp() with a single neighbouring cell returns that cell
  mask <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0,
                      ymax = 1, crs = "EPSG:4326")
  terra::values(mask) <- 1:2
  for (seed in 1:10) {
    set.seed(seed)
    testthat::expect_equal(
      tenm::cells2samp(data = data.frame(lon = 0.5, lat = 0.5),
                       longitude = "lon", latitude = "lat", buffer_ngbs = 1,
                       raster_mask = mask, n_bg = 1, progress = FALSE), 2)
  }
})

testthat::test_that("#13: pROC() honours rseed and keeps the caller's RNG state", {
  set.seed(1)
  bgv <- runif(5000)
  tv <- runif(30, .3, 1)
  # rseed = FALSE: the caller's stream continues (it is not reset)
  set.seed(111)
  r111 <- tenm::pROC(bgv, tv, n_iter = 20, rseed = FALSE, sub_sample = FALSE)
  a <- runif(1)
  set.seed(999)
  r999 <- tenm::pROC(bgv, tv, n_iter = 20, rseed = FALSE, sub_sample = FALSE)
  b <- runif(1)
  testthat::expect_false(identical(a, b))
  testthat::expect_false(identical(r111$pROC_results, r999$pROC_results))
  # rseed = TRUE: reproducible, and .Random.seed is left unchanged
  set.seed(222)
  before <- .Random.seed
  t1 <- tenm::pROC(bgv, tv, n_iter = 20, rseed = TRUE)
  testthat::expect_identical(.Random.seed, before)
  set.seed(222)
  t2 <- tenm::pROC(bgv, tv, n_iter = 20, rseed = TRUE)
  testthat::expect_identical(t1$pROC_results, t2$pROC_results)
})

testthat::test_that("#14: ellipsoid_omr() ignores a global object called in_etest", {
  set.seed(14)
  env_train <- data.frame(bio_01 = rnorm(60, 20, 3),
                          bio_12 = rnorm(60, 1000, 150))
  env_bg <- data.frame(bio_01 = runif(2000, 5, 40),
                       bio_12 = runif(2000, 200, 2500))
  ref <- tenm::ellipsoid_omr(env_data = env_train, env_test = NULL,
                             env_bg = env_bg, cf_level = 0.975, mve = FALSE)
  assign("in_etest", "anything", envir = globalenv())
  on.exit(rm("in_etest", envir = globalenv()), add = TRUE)
  got <- tenm::ellipsoid_omr(env_data = env_train, env_test = NULL,
                             env_bg = env_bg, cf_level = 0.975, mve = FALSE)
  testthat::expect_equal(got, ref)
  testthat::expect_false("pval_bin" %in% names(got))
})

testthat::test_that("#15: tenm_selection() passes its pROC settings to pROC()", {
  set.seed(15)
  tdf <- data.frame(bio_01 = rnorm(80, 20, 3), bio_12 = rnorm(80, 1000, 150),
                    trian_test = rep(c("Train", "Test"), c(60, 20)))
  sp_bg <- list(temporal_df = tdf, layers_ext = ".tif$",
                env_bg = data.frame(bio_01 = runif(3000, 5, 40),
                                    bio_12 = runif(3000, 200, 2500)))
  class(sp_bg) <- "sp.temporal.bg"
  utils::capture.output(
    sel <- tenm::tenm_selection(sp_bg, omr_criteria = 0.1,
                                ellipsoid_level = 0.975,
                                vars2fit = c("bio_01", "bio_12"),
                                nvars_to_fit = 2, mve = FALSE, proc = TRUE,
                                RandomPercent = 30, sub_sample = FALSE,
                                sub_sample_size = 200, NoOfIteration = 9,
                                E_percent = 10, parallel = FALSE))
  train <- tdf[tdf$trian_test == "Train", 1:2]
  test <- tdf[tdf$trian_test == "Test", 1:2]
  direct <- tenm::ellipsoid_omr(train, test, sp_bg$env_bg, cf_level = 0.975,
                                mve = FALSE, proc = TRUE, proc_iter = 9,
                                rseed = TRUE, E_percent = 10,
                                boost_percent = 30, sub_sample = FALSE,
                                sub_sample_size = 200)
  default <- tenm::ellipsoid_omr(train, test, sp_bg$env_bg, cf_level = 0.975,
                                 mve = FALSE, proc = TRUE, proc_iter = 9,
                                 rseed = TRUE, sub_sample = FALSE)
  cols <- c("pval_proc", "env_bg_paucratio", "env_bg_auc")
  testthat::expect_equal(unlist(sel$mods_table[1, cols]), unlist(direct[1, cols]))
  testthat::expect_false(isTRUE(all.equal(unlist(direct[1, cols]),
                                          unlist(default[1, cols]))))
})

testthat::test_that("#16: predict() accepts plot = FALSE, names bad variables and defaults layers_ext", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  dir70 <- system.file("extdata/bio_1970_2000", package = "tenm")
  bio_01 <- terra::rast(file.path(dir70, "bio_01.tif"))
  bio_12 <- terra::rast(file.path(dir70, "bio_12.tif"))
  set.seed(16)
  cells <- which(!is.na(terra::values(bio_01)[, 1]) &
                   !is.na(terra::values(bio_12)[, 1]))
  cells <- sample(cells, 40)
  train <- data.frame(bio_01 = terra::values(bio_01)[cells, 1],
                      bio_12 = terra::values(bio_12)[cells, 1],
                      trian_test = "Train")
  mod_sel <- list(temporal_df = train, layers_ext = ".tif$",
                  env_bg = train[, c("bio_01", "bio_12")],
                  mods_table = data.frame(fitted_vars = "bio_01,bio_12"))
  class(mod_sel) <- "sp.temporal.selection"
  utils::capture.output({
    p1 <- predict(mod_sel, model_variables = c("bio_01", "bio_12"),
                  layers_path = dir70, layers_ext = ".tif$", mve = FALSE,
                  plot = FALSE)
    p2 <- predict(mod_sel, model_variables = c("bio_01", "bio_12"),
                  layers_path = dir70, mve = FALSE)
  })
  testthat::expect_s4_class(p1, "SpatRaster")
  testthat::expect_equal(terra::values(p2), terra::values(p1),
                         ignore_attr = TRUE)
  testthat::expect_error(
    predict(mod_sel, model_variables = c("bio_01", "bio_99"),
            layers_path = dir70, layers_ext = ".tif$"),
    "Not valid variable names: bio_99")
})

testthat::test_that("#20: non_pred_test_ids is NA when no test point is omitted", {
  set.seed(20)
  env_train <- data.frame(bio_01 = rnorm(60, 20, 3),
                          bio_12 = rnorm(60, 1000, 150))
  env_test <- data.frame(bio_01 = c(20, 20.5, 19.5, 21, 19),
                         bio_12 = c(1000, 990, 1010, 1020, 980))
  env_bg <- data.frame(bio_01 = runif(2000, 5, 40),
                       bio_12 = runif(2000, 200, 2500))
  r_in <- tenm::ellipsoid_omr(env_train, env_test, env_bg, cf_level = 0.975,
                              mve = FALSE)
  testthat::expect_equal(r_in$om_rate_test, 0)
  testthat::expect_true(is.na(r_in$non_pred_test_ids))
  r_out <- tenm::ellipsoid_omr(env_train,
                               rbind(env_test, data.frame(bio_01 = 60,
                                                          bio_12 = 5000)),
                               env_bg, cf_level = 0.975, mve = FALSE)
  testthat::expect_equal(as.character(r_out$non_pred_test_ids), "6")
})

testthat::test_that("#22: pROC() handles one unique test value and integer suitabilities", {
  set.seed(1)
  bgv <- runif(5000)
  a <- tenm::pROC(bgv, rep(0.5, 10), n_iter = 5, rseed = TRUE)$pROC_summary
  testthat::expect_true(all(is.finite(a)))
  b <- tenm::pROC(0:100, c(40, 60, 80, 90), n_iter = 5, rseed = TRUE,
                  sub_sample = FALSE)$pROC_summary
  testthat::expect_true(all(is.finite(b)))
  set.seed(22)
  env_train <- data.frame(bio_01 = rnorm(60, 20, 3),
                          bio_12 = rnorm(60, 1000, 150))
  env_test <- data.frame(bio_01 = rep(20, 5), bio_12 = rep(1000, 5))
  env_bg <- data.frame(bio_01 = runif(2000, 5, 40),
                       bio_12 = runif(2000, 200, 2500))
  r <- tenm::ellipsoid_omr(env_train, env_test, env_bg, cf_level = 0.975,
                           mve = FALSE, proc = TRUE, proc_iter = 5)
  testthat::expect_false(anyNA(r[, c("pval_proc", "env_bg_paucratio",
                                     "env_bg_auc")]))
})

testthat::test_that("#24: predict() refits with the mve and level used in selection", {
  dir70 <- system.file("extdata/bio_1970_2000", package = "tenm")
  lyr <- terra::rast(file.path(dir70, c("bio_01.tif", "bio_12.tif")))
  names(lyr) <- c("bio_01", "bio_12")
  v <- terra::values(lyr)
  ok_rows <- stats::complete.cases(v)
  set.seed(24)
  cells <- sample(which(ok_rows), 80)
  tdf <- data.frame(bio_01 = v[cells, 1], bio_12 = v[cells, 2],
                    trian_test = rep(c("Train", "Test"), c(60, 20)))
  tdf[1:6, "bio_12"] <- tdf[1:6, "bio_12"] * 3
  sp_bg <- list(temporal_df = tdf, layers_ext = ".tif$",
                env_bg = as.data.frame(v[ok_rows, ]))
  class(sp_bg) <- "sp.temporal.bg"
  utils::capture.output(
    sel <- tenm::tenm_selection(sp_bg, omr_criteria = 0.1,
                                ellipsoid_level = 0.95,
                                vars2fit = c("bio_01", "bio_12"),
                                nvars_to_fit = 2, mve = FALSE, proc = FALSE,
                                parallel = FALSE))
  testthat::expect_false(sel$mve)
  testthat::expect_equal(sel$ellipsoid_level, 0.95)
  train <- tdf[tdf$trian_test == "Train", c("bio_01", "bio_12")]
  ref_mod <- tenm::cov_center(train, mve = FALSE, level = 0.95,
                              vars = c("bio_01", "bio_12"))
  ref <- rep(NA_real_, nrow(v))
  ref[ok_rows] <- exp(-0.5 * stats::mahalanobis(v[ok_rows, ], ref_mod$centroid,
                                                ref_mod$covariance))
  utils::capture.output(p <- predict(sel, layers = lyr), type = "message")
  testthat::expect_equal(terra::values(p)[, 1], ref)
  # An explicit mve still overrides the stored setting
  set.seed(1)
  utils::capture.output(p_mve <- predict(sel, layers = lyr, mve = TRUE),
                        type = "message")
  testthat::expect_false(isTRUE(all.equal(terra::values(p_mve)[, 1], ref)))
})

testthat::test_that("#26: missing test omission stays NA, failed fits are reported, columns are stable", {
  set.seed(26)
  env_train <- data.frame(bio_01 = rnorm(60, 20, 3),
                          bio_12 = rnorm(60, 1000, 150),
                          bio_03 = rnorm(60, 50, 5))
  env_test <- data.frame(bio_01 = rnorm(20, 20, 3),
                         bio_12 = rnorm(20, 1000, 150), bio_03 = NA_real_)
  env_bg <- data.frame(bio_01 = runif(2000, 5, 40),
                       bio_12 = runif(2000, 200, 2500),
                       bio_03 = runif(2000, 20, 80))
  env_train$bio_99 <- 2 * env_train$bio_01
  env_test$bio_99 <- 2 * env_test$bio_01
  env_bg$bio_99 <- 2 * env_bg$bio_01
  run <- function(crit) {
    out <- NULL
    utils::capture.output(
      out <- tenm::ellipsoid_selection(env_train, env_test,
                                       env_vars = c("bio_01", "bio_12",
                                                    "bio_03", "bio_99"),
                                       nvarstest = 2, level = 0.975,
                                       mve = FALSE, env_bg = env_bg,
                                       omr_criteria = crit, parallel = FALSE,
                                       proc = TRUE, proc_iter = 5),
      type = "message")
    out
  }
  utils::capture.output({
    sel <- run(0.1)
    sel_none <- run(-1)
  })
  no_test <- grepl("bio_03", sel$fitted_vars)
  ok <- sel$status == "ok"
  testthat::expect_true(all(is.na(sel$om_rate_test[no_test])))
  testthat::expect_true(max(sel$rank_by_omr_train_test[ok & !no_test]) <
                          min(sel$rank_by_omr_train_test[no_test]))
  failed <- sel[sel$status == "fit_failed", ]
  testthat::expect_equal(failed$fitted_vars, "bio_01,bio_99")
  testthat::expect_match(failed$fit_error, "singular")
  testthat::expect_equal(nrow(sel), 6)
  testthat::expect_identical(names(sel), names(sel_none))
  testthat::expect_true("rank_omr_aucratio" %in% names(sel_none))
  # ellipsoid_omr() reports a failed fit instead of returning NULL
  r <- tenm::ellipsoid_omr(env_train[, c("bio_01", "bio_99")], NULL,
                           env_bg[, c("bio_01", "bio_99")],
                           cf_level = 0.975, mve = FALSE)
  testthat::expect_equal(r$status, "fit_failed")
})

testthat::test_that("#27: bg_prevalence is computed over the background only", {
  set.seed(27)
  env_train <- data.frame(bio_01 = rnorm(60, 20, 3),
                          bio_12 = rnorm(60, 1000, 150))
  env_test <- data.frame(bio_01 = rnorm(20, 20, 3),
                         bio_12 = rnorm(20, 1000, 150))
  env_bg <- data.frame(bio_01 = runif(300, 5, 40),
                       bio_12 = runif(300, 200, 2500))
  r <- tenm::ellipsoid_omr(env_train, env_test, env_bg, cf_level = 0.975,
                           mve = FALSE)
  m <- tenm::cov_center(env_train, mve = FALSE, level = 0.975, vars = 1:2)
  in_bg <- tenm::inEllipsoid(m$centroid, m$covariance, env_bg,
                             level = 0.975)$in_Ellipsoid
  in_te <- tenm::inEllipsoid(m$centroid, m$covariance, env_test,
                             level = 0.975)$in_Ellipsoid
  testthat::expect_equal(r$bg_prevalence, mean(in_bg))
  testthat::expect_equal(r$pval_bin,
                         stats::binom.test(sum(in_te), length(in_te),
                                           p = mean(in_bg),
                                           alternative = "greater")$p.value)
})

testthat::test_that("#9: clean_dup() stops when by_mask = TRUE and raster_mask is not a SpatRaster", {
  lyr <- make_year_layers("2000")
  mask <- terra::rast(file.path(lyr$root, "2000", "bio_01.tif"))
  occ <- data.frame(lon = c(1.5, 1.6, 5.5), lat = c(1.5, 1.6, 5.5), id = 1:3)
  testthat::expect_error(
    tenm::clean_dup(occ, "lon", "lat", by_mask = TRUE, raster_mask = mask[1],
                    progress = FALSE), "SpatRaster")
  testthat::expect_error(
    tenm::clean_dup(occ, "lon", "lat", by_mask = TRUE, raster_mask = NULL,
                    progress = FALSE), "SpatRaster")
  testthat::expect_equal(
    nrow(tenm::clean_dup(occ, "lon", "lat", by_mask = TRUE, raster_mask = mask,
                         progress = FALSE)), 2)
})
