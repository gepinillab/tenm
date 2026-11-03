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
