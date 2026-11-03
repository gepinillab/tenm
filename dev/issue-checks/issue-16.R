# Issue #16: predict(): plot = FALSE errors, the bad-variable message prints
# NULL, and the default layers_ext = NULL errors.
# Uses the bundled 1970-2000 layers and a small synthetic training set.
# Run with: Rscript dev/issue-checks/issue-16.R
cat("Issue #16: predict(): plot = FALSE errors, the bad-variable message prints NULL, and the default layers_ext = NULL errors\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")
grDevices::pdf(NULL)  # predict() may draw a plot; discard it

dir70 <- system.file("extdata/bio_1970_2000", package = "tenm")
bio_01 <- terra::rast(file.path(dir70, "bio_01.tif"))
bio_12 <- terra::rast(file.path(dir70, "bio_12.tif"))
set.seed(16)
cells <- which(!is.na(terra::values(bio_01)[, 1]) & !is.na(terra::values(bio_12)[, 1]))
cells <- sample(cells, 40)
train <- data.frame(bio_01 = terra::values(bio_01)[cells, 1],
                    bio_12 = terra::values(bio_12)[cells, 1],
                    trian_test = "Train")
mod_sel <- list(temporal_df = train, layers_ext = ".tif$",
                env_bg = train[, c("bio_01", "bio_12")],
                mods_table = data.frame(fitted_vars = "bio_01,bio_12"))
class(mod_sel) <- "sp.temporal.selection"

# Record the plot argument ellipsoid_projection() receives.
rec <- new.env()
rec$plot <- c()
invisible(suppressMessages(trace("ellipsoid_projection",
  tracer = quote(assign("plot", c(get("plot", envir = rec), plot), envir = rec)),
  where = asNamespace("tenm"), print = FALSE)))

quiet <- function(expr) {
  out <- NULL
  invisible(utils::capture.output(out <- tryCatch(expr, error = function(e) e)))
  out
}
# A) plot = FALSE
p_a <- quiet(predict(mod_sel, model_variables = c("bio_01", "bio_12"), layers_path = dir70,
                     layers_ext = ".tif$", mve = FALSE, plot = FALSE))
plot_a <- rec$plot
# B) an invalid variable name
p_b <- quiet(predict(mod_sel, model_variables = c("bio_01", "bio_99"), layers_path = dir70,
                     layers_ext = ".tif$", mve = FALSE))
# C) layers_ext left at its default
p_c <- quiet(predict(mod_sel, model_variables = c("bio_01", "bio_12"), layers_path = dir70,
                     mve = FALSE))
suppressMessages(untrace("ellipsoid_projection", where = asNamespace("tenm")))

ok_a <- inherits(p_a, "SpatRaster") && identical(plot_a, FALSE)
ok_b <- inherits(p_b, "error") && grepl("bio_99", conditionMessage(p_b)) &&
  !grepl("NULL", conditionMessage(p_b))
ok_c <- inherits(p_c, "SpatRaster") && inherits(p_a, "SpatRaster") &&
  isTRUE(all.equal(terra::values(p_c), terra::values(p_a), check.attributes = FALSE))
msg <- function(x) if (inherits(x, "error")) paste0("error \"", gsub("\n", " ", conditionMessage(x)), "\"") else
  class(x)[1]

cat("EXPECTED: A) predict(plot = FALSE) returns a SpatRaster and no plot is drawn;",
    "B) model_variables = c(\"bio_01\", \"bio_99\") errors naming bio_99;",
    "C) predict() without layers_ext uses the object's layers_ext and returns the same raster as A\n")
cat("GOT:      A)", msg(p_a), if (!inherits(p_a, "error")) paste("| plot =", plot_a), "\n")
cat("          B)", msg(p_b), "\n")
cat("          C)", msg(p_c), "\n")
cat("RESULT:", if (ok_a && ok_b && ok_c) "PASS" else "FAIL", "\n")
