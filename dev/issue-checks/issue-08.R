# Issue #8: clean_dup_by_date() fails under future multisession plans
# (globals whitelist).
# Synthetic year folders 2000 and 2001; threshold and raster_mask are passed
# as variables, as in a script that sets plan("multisession") early.
# Run with: Rscript dev/issue-checks/issue-08.R
cat("Issue #8: clean_dup_by_date() fails under future multisession plans (globals whitelist)\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

root <- file.path(tempdir(), "issue08_layers")
unlink(root, recursive = TRUE)
dir.create(root)
for (f in c("2000", "2001")) {
  dir.create(file.path(root, f))
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                   crs = "EPSG:4326", vals = 1:100)
  names(r) <- "bio_01"
  terra::writeRaster(r, file.path(root, f, "bio_01.tif"))
}
occ <- data.frame(lon = c(1.5, 1.6, 5.5, 1.5, 1.6, 8.5), lat = c(1.5, 1.6, 5.5, 1.5, 1.6, 8.5),
                  year = rep(c(2000, 2001), each = 3))
st <- sp_temporal_data(occ, "lon", "lat", "year", "y", "y", root)
tenm_mask <- terra::rast(file.path(root, "2000", "bio_01.tif"))
th <- 0.5

run <- function(expr) tryCatch(nrow(expr$temporal_df), error = function(e) {
  paste("error:", gsub("\n", " ", conditionMessage(e)))
})
future::plan("sequential")
seq_dist <- run(clean_dup_by_date(st, threshold = th, progress = FALSE))
seq_mask <- run(clean_dup_by_date(st, threshold = 1, by_mask = TRUE,
                                  raster_mask = tenm_mask, progress = FALSE))
future::plan("multisession", workers = 2)
ms_dist <- run(clean_dup_by_date(st, threshold = th, progress = FALSE))
ms_mask <- run(clean_dup_by_date(st, threshold = 1, by_mask = TRUE,
                                 raster_mask = tenm_mask, progress = FALSE))
future::plan("sequential")

ok <- identical(ms_dist, seq_dist) && identical(ms_mask, seq_mask) &&
  is.numeric(seq_dist) && is.numeric(seq_mask)
cat("EXPECTED: under plan(multisession), threshold = th gives", seq_dist,
    "rows and raster_mask = tenm_mask gives", seq_mask, "rows, as under plan(sequential)\n")
cat("GOT:      threshold = th:", ms_dist, "\n")
cat("          raster_mask = tenm_mask:", ms_mask, "\n")
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
