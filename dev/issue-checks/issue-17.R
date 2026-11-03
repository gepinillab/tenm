# Issue #17: sp_temporal_data(): tibble input errors, and the documented
# "mdy" format is rejected.
# Synthetic year folders 2000 and 2001.
# Run with: Rscript dev/issue-checks/issue-17.R
cat("Issue #17: sp_temporal_data(): tibble input errors, and the documented \"mdy\" format is rejected\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

root <- file.path(tempdir(), "issue17_layers")
unlink(root, recursive = TRUE)
dir.create(root)
for (f in c("2000", "2001")) {
  dir.create(file.path(root, f))
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                   crs = "EPSG:4326", vals = 1:100)
  names(r) <- "bio_01"
  terra::writeRaster(r, file.path(root, f, "bio_01.tif"))
}
occ <- data.frame(lon = c(1.5, 5.5, 8.5), lat = c(1.5, 5.5, 8.5), year = c(2000, 2001, 2001),
                  date = c("06-15-2000", "01-20-2001", "11-30-2001"))

run <- function(expr) tryCatch(suppressMessages(expr), error = function(e) {
  paste("error:", gsub("\n", " ", conditionMessage(e)))
})
tb <- run(sp_temporal_data(tibble::as_tibble(occ), "lon", "lat", "year", "y", "y", root))
df <- sp_temporal_data(occ, "lon", "lat", "year", "y", "y", root)
md <- run(sp_temporal_data(occ, "lon", "lat", "date", "mdy", "y", root))

ok_tb <- is.list(tb) && isTRUE(all.equal(as.data.frame(tb$temporal_df), df$temporal_df,
                                         check.attributes = FALSE))
ok_md <- is.list(md) && identical(basename(md$temporal_df$layers_path), c("2000", "2001", "2001"))
cat("EXPECTED: tibble input gives the same temporal_df as a data.frame;",
    "occ_date_format = 'mdy' reads 06-15-2000 etc. and matches folders 2000, 2001, 2001\n")
cat("GOT:      tibble:", if (is.list(tb)) paste(nrow(tb$temporal_df), "rows",
                                                  if (ok_tb) "(same as data.frame)" else "(differs)") else tb, "\n")
cat("          mdy:   ", if (is.list(md)) paste("folders", paste(basename(md$temporal_df$layers_path),
                                                                collapse = ", ")) else md, "\n")
cat("RESULT:", if (ok_tb && ok_md) "PASS" else "FAIL", "\n")
