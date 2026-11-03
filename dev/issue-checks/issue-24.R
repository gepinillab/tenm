# Issue #24: predict() refits the model with its own mve/level instead of the
# selected model's settings.
# Uses the bundled 1970-2000 layers and a small synthetic training set.
# Run with: Rscript dev/issue-checks/issue-24.R
cat("Issue #24: predict() refits the model with its own mve/level instead of the selected model's settings\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")
grDevices::pdf(NULL)  # predict() may draw a plot; discard it

dir70 <- system.file("extdata/bio_1970_2000", package = "tenm")
lyr <- terra::rast(file.path(dir70, c("bio_01.tif", "bio_12.tif")))
names(lyr) <- c("bio_01", "bio_12")
v <- terra::values(lyr)
set.seed(24)
cells <- sample(which(stats::complete.cases(v)), 80)
# Training points with a few outliers, so mve = TRUE and mve = FALSE differ.
tdf <- data.frame(bio_01 = v[cells, 1], bio_12 = v[cells, 2],
                  trian_test = rep(c("Train", "Test"), c(60, 20)))
tdf[1:6, "bio_12"] <- tdf[1:6, "bio_12"] * 3
sp_bg <- list(temporal_df = tdf, sp_date_var = "year", lon_lat_vars = c("x", "y"),
              layers_ext = ".tif$", env_bg = as.data.frame(v[stats::complete.cases(v), ]))
class(sp_bg) <- "sp.temporal.bg"

invisible(utils::capture.output(
  sel <- tenm_selection(sp_bg, omr_criteria = 0.1, ellipsoid_level = 0.95,
                        vars2fit = c("bio_01", "bio_12"), nvars_to_fit = 2,
                        mve = FALSE, proc = FALSE, parallel = FALSE)))

# Record the mve and level that cov_center() receives inside predict().
rec <- new.env()
rec$args <- list()
invisible(suppressMessages(trace("cov_center",
  tracer = quote(assign("args", c(get("args", envir = rec), list(list(mve = mve, level = level))),
                        envir = rec)),
  where = asNamespace("tenm"), print = FALSE)))
invisible(utils::capture.output(p <- predict(sel, layers = lyr)))
suppressMessages(untrace("cov_center", where = asNamespace("tenm")))
got <- rec$args[[1]]

# Reference: projection of the mve = FALSE ellipsoid that selection evaluated.
train <- tdf[tdf$trian_test == "Train", c("bio_01", "bio_12")]
ref_mod <- cov_center(train, mve = FALSE, level = 0.95, vars = c("bio_01", "bio_12"))
ok_rows <- stats::complete.cases(v)
ref <- rep(NA_real_, nrow(v))
ref[ok_rows] <- exp(-0.5 * stats::mahalanobis(v[ok_rows, ], ref_mod$centroid, ref_mod$covariance))
max_diff <- max(abs(terra::values(p)[, 1] - ref), na.rm = TRUE)

ok <- isFALSE(got$mve) && isTRUE(all.equal(got$level, 0.95)) && max_diff < 1e-10

cat("EXPECTED: after tenm_selection(mve = FALSE, ellipsoid_level = 0.95), predict() without",
    "mve/level refits with mve = FALSE, level = 0.95 and reproduces the evaluated ellipsoid",
    "(max |difference| < 1e-10)\n")
cat("GOT:      predict() refit with mve =", got$mve, "| level =", got$level,
    "| max |difference| from the evaluated ellipsoid =", format(max_diff, digits = 4), "\n")
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
