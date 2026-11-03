# Issue #3: predict(): layers are stacked in file order, not model-variable
# order, giving wrong projections.
# Uses the bundled 1970-2000 layers and a small synthetic training set.
# Run with: Rscript dev/issue-checks/issue-03.R
cat("Issue #3: predict(): layers are stacked in file order, not model-variable order, giving wrong projections\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")
grDevices::pdf(NULL)  # predict() draws a plot; discard it

dir70 <- system.file("extdata/bio_1970_2000", package = "tenm")
bio_01 <- terra::rast(file.path(dir70, "bio_01.tif"))
bio_12 <- terra::rast(file.path(dir70, "bio_12.tif"))
names(bio_01) <- "bio_01"
names(bio_12) <- "bio_12"

# Synthetic "sp.temporal.selection" object: 40 training points taken from
# non-NA cells of the 1970-2000 layers.
set.seed(3)
cells <- which(!is.na(terra::values(bio_01)[, 1]) & !is.na(terra::values(bio_12)[, 1]))
cells <- sample(cells, 40)
train <- data.frame(bio_01 = terra::values(bio_01)[cells, 1],
                    bio_12 = terra::values(bio_12)[cells, 1],
                    trian_test = "Train")
mod_sel <- list(temporal_df = train,
                env_bg = train[, c("bio_01", "bio_12")],
                mods_table = data.frame(fitted_vars = "bio_01,bio_12"))
class(mod_sel) <- "sp.temporal.selection"

# Reference suitability, computed by hand with the layers in model order.
mod <- cov_center(train[, c("bio_01", "bio_12")], mve = FALSE, level = 0.975,
                  vars = c("bio_01", "bio_12"))
vals <- cbind(bio_01 = terra::values(bio_01)[, 1], bio_12 = terra::values(bio_12)[, 1])
ref <- rep(NA_real_, nrow(vals))
ok_rows <- stats::complete.cases(vals)
ref[ok_rows] <- exp(-0.5 * stats::mahalanobis(vals[ok_rows, ], center = mod$centroid,
                                              cov = mod$covariance))

max_diff <- function(p) {
  v <- terra::values(p)[, 1]
  max(abs(v - ref), na.rm = TRUE)
}
safe <- function(expr) tryCatch(expr, error = function(e) {
  cat("  error:", conditionMessage(e), "\n")
  NULL
})

invisible(utils::capture.output(suppressMessages({
  # layers_path branch: same model, variables listed in two orders
  p_path_12 <- safe(predict(mod_sel, model_variables = c("bio_01", "bio_12"),
                            layers_path = dir70, layers_ext = ".tif$", mve = FALSE))
  p_path_21 <- safe(predict(mod_sel, model_variables = c("bio_12", "bio_01"),
                            layers_path = dir70, layers_ext = ".tif$", mve = FALSE))
  # SpatRaster branch: layers stacked in a different order than the model
  p_rast <- safe(predict(mod_sel, model_variables = c("bio_01", "bio_12"),
                         layers = c(bio_12, bio_01), mve = FALSE))
  # list branch: same, inside a list
  p_list <- safe(predict(mod_sel, model_variables = c("bio_01", "bio_12"),
                         layers = list(c(bio_12, bio_01)), mve = FALSE))
}), type = "output"))
if (is.list(p_list)) p_list <- p_list[[1]]

diffs <- c(
  path_bio01_bio12 = if (is.null(p_path_12)) NA else max_diff(p_path_12),
  path_bio12_bio01 = if (is.null(p_path_21)) NA else max_diff(p_path_21),
  spatraster_reordered = if (is.null(p_rast)) NA else max_diff(p_rast),
  list_reordered = if (is.null(p_list)) NA else max_diff(p_list))

cat("\nEXPECTED: every projection equals the reference (layers matched by name);",
    "max |difference| < 1e-10 for all four calls\n")
cat("GOT:      max |difference| =",
    paste(names(diffs), format(diffs, digits = 4), sep = ": ", collapse = " | "), "\n")
cat("RESULT:", if (all(!is.na(diffs) & diffs < 1e-10)) "PASS" else "FAIL", "\n")
