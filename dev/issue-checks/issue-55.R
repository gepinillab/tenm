# Issue #55: predict() hard-codes size = 2, so predict(..., size = 3) fails
# with "formal argument "size" matched by multiple actual arguments".
# Synthetic object and rasters; plots are drawn to a temporary PNG.
# Run with: Rscript dev/issue-checks/issue-55.R
cat("Issue #55: predict() hard-codes size = 2\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

set.seed(55)
object <- structure(
  list(temporal_df = data.frame(x = rnorm(30), y = rnorm(30), trian_test = "Train"),
       env_bg = data.frame(x = 0, y = 0),
       mods_table = data.frame(fitted_vars = "x,y")),
  class = "sp.temporal.selection")
layers <- terra::rast(nrows = 2, ncols = 2, nlyrs = 2,
                      xmin = 0, xmax = 2, ymin = 0, ymax = 2)
names(layers) <- c("x", "y")
terra::values(layers) <- matrix(rnorm(8), ncol = 2)
layer_dir <- tempfile()
dir.create(layer_dir)
terra::writeRaster(layers[[1]], file.path(layer_dir, "x.tif"))
terra::writeRaster(layers[[2]], file.path(layer_dir, "y.tif"))

try_predict <- function(...) {
  grDevices::png(tempfile(fileext = ".png"))
  on.exit(grDevices::dev.off())
  tryCatch({
    out <- NULL
    utils::capture.output(suppressMessages(out <- predict(object, mve = FALSE, ...)))
    if (is.list(out)) "ok" else if (methods::is(out, "SpatRaster")) "ok" else "unexpected output"
  }, error = function(e) paste("error:", conditionMessage(e)))
}
res <- c(
  raster = try_predict(layers = layers, size = 3),
  raster_plot = try_predict(layers = layers, plot = TRUE, size = 3),
  list = try_predict(layers = list(layers), size = 3),
  dir = try_predict(layers_path = layer_dir, layers_ext = "\\.tif$", size = 3),
  default = try_predict(layers = layers)  # existing call, no size
)
cat("EXPECTED: predict(..., size = 3) works for a raster, a raster list and a\n",
    "         directory; calls without size still work\n")
cat("GOT:\n")
print(res)
cat("\nRESULT:", if (all(res == "ok")) "PASS" else "FAIL", "\n")
