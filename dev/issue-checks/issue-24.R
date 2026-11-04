# Issue #24: predict() refits the model with its own mve/level instead of the
# selected model's settings.
# Uses the bundled 1970-2000 layers and a small synthetic training set.
# Checks that tenm_selection() stores the fitted centroid and covariance of
# each model that passes selection, that they equal the ellipsoid evaluated
# during selection, that predict() uses them (repeated calls with mve = TRUE
# are identical), that an mve = FALSE selection is reproduced, and that an
# object without stored ellipsoids still predicts (refit fallback).
# Run with: Rscript dev/issue-checks/issue-24.R
cat("Issue #24: predict() refits the model with its own mve/level instead of the selected model's settings\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")
grDevices::pdf(NULL)  # predict() may draw a plot; discard it

dir70 <- system.file("extdata/bio_1970_2000", package = "tenm")
lyr <- terra::rast(file.path(dir70, c("bio_01.tif", "bio_12.tif")))
names(lyr) <- c("bio_01", "bio_12")
v <- terra::values(lyr)
ok_rows <- stats::complete.cases(v)
set.seed(24)
cells <- sample(which(ok_rows), 80)
# Training points with a few outliers, so mve = TRUE and mve = FALSE differ.
tdf <- data.frame(bio_01 = v[cells, 1], bio_12 = v[cells, 2],
                  trian_test = rep(c("Train", "Test"), c(60, 20)))
tdf[1:6, "bio_12"] <- tdf[1:6, "bio_12"] * 3
sp_bg <- list(temporal_df = tdf, sp_date_var = "year", lon_lat_vars = c("x", "y"),
              layers_ext = ".tif$", env_bg = as.data.frame(v[ok_rows, ]))
class(sp_bg) <- "sp.temporal.bg"
vars <- c("bio_01", "bio_12")

# (a) mve = TRUE, with a diffuse synthetic training cloud so that the
# stochastic MASS::cov.mve fit differs between refits.
rng <- apply(v[ok_rows, ], 2, range)
tdf_a <- data.frame(bio_01 = stats::runif(200, rng[1, 1], rng[2, 1]),
                    bio_12 = stats::runif(200, rng[1, 2], rng[2, 2]),
                    trian_test = rep(c("Train", "Test"), c(150, 50)))
sp_bg_a <- sp_bg
sp_bg_a$temporal_df <- tdf_a
# Record the ellipsoid that cov_center() returns during
# selection (the evaluated one), then compare it with the stored one.
rec <- new.env()
rec$fits <- list()
invisible(suppressMessages(trace("cov_center",
  exit = quote(assign("fits", c(get("fits", envir = rec), list(returnValue())), envir = rec)),
  where = asNamespace("tenm"), print = FALSE)))
invisible(utils::capture.output(
  sel <- tenm_selection(sp_bg_a, omr_criteria = 0.2, ellipsoid_level = 0.975,
                        vars2fit = vars, nvars_to_fit = 2,
                        mve = TRUE, proc = FALSE, parallel = FALSE)))
suppressMessages(untrace("cov_center", where = asNamespace("tenm")))
evaluated <- rec$fits[[1]]
stored <- if (is.list(sel$ellipsoids)) sel$ellipsoids[["bio_01,bio_12"]] else NULL
same_ell <- !is.null(stored) && isTRUE(all.equal(stored$centroid, evaluated$centroid)) &&
  isTRUE(all.equal(stored$covariance, evaluated$covariance))

# Repeated predictions with mve = TRUE, with the RNG in different states
set.seed(1)
invisible(utils::capture.output(p1 <- predict(sel, layers = lyr), type = "message"))
set.seed(2)
invisible(utils::capture.output(p2 <- predict(sel, layers = lyr), type = "message"))
repeat_same <- identical(terra::values(p1), terra::values(p2))
ref_a <- rep(NA_real_, nrow(v))
ref_a[ok_rows] <- exp(-0.5 * stats::mahalanobis(v[ok_rows, ], evaluated$centroid,
                                                evaluated$covariance))
diff_a <- max(abs(terra::values(p1)[, 1] - ref_a), na.rm = TRUE)

# (b) mve = FALSE, ellipsoid_level = 0.95: predict() reproduces that ellipsoid
invisible(utils::capture.output(
  sel_f <- tenm_selection(sp_bg, omr_criteria = 0.2, ellipsoid_level = 0.95,
                          vars2fit = vars, nvars_to_fit = 2,
                          mve = FALSE, proc = FALSE, parallel = FALSE)))
train <- tdf[tdf$trian_test == "Train", vars]
ref_mod <- cov_center(train, mve = FALSE, level = 0.95, vars = vars)
ref_b <- rep(NA_real_, nrow(v))
ref_b[ok_rows] <- exp(-0.5 * stats::mahalanobis(v[ok_rows, ], ref_mod$centroid,
                                                ref_mod$covariance))
invisible(utils::capture.output(pf <- predict(sel_f, layers = lyr), type = "message"))
diff_b <- max(abs(terra::values(pf)[, 1] - ref_b), na.rm = TRUE)

# (c) an object saved before the ellipsoids were stored: refit fallback
old <- sel_f
old$ellipsoids <- NULL
pc <- tryCatch({
  invisible(utils::capture.output(x <- predict(old, layers = lyr), type = "message"))
  x
}, error = function(e) NULL)
diff_c <- if (is.null(pc)) NA else max(abs(terra::values(pc)[, 1] - ref_b), na.rm = TRUE)

ok <- same_ell && repeat_same && diff_a < 1e-10 && diff_b < 1e-10 &&
  !is.na(diff_c) && diff_c < 1e-10

cat("EXPECTED: (a) mve = TRUE: the object stores the centroid and covariance of the passing",
    "model, equal to the ellipsoid evaluated in selection; two predict() calls are identical",
    "and project that ellipsoid; (b) mve = FALSE, level = 0.95 is reproduced (< 1e-10);",
    "(c) an object without stored ellipsoids still predicts by refitting (< 1e-10)\n")
cat("GOT:      (a) stored ellipsoid:", if (is.null(stored)) "none" else
  paste("present, equal to evaluated =", same_ell),
  "| repeated predict() identical =", repeat_same,
  "| max |difference| from evaluated =", format(diff_a, digits = 4), "\n")
cat("          (b) max |difference| =", format(diff_b, digits = 4), "\n")
cat("          (c)", if (is.null(pc)) "error" else
  paste("max |difference| =", format(diff_c, digits = 4)), "\n")
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
