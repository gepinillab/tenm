# Issue #7: ex_by_date()/bg_by_date(): first year's cell IDs are reused for
# all years, giving wrong background values when grids differ.
# Synthetic layers: the 1940 layers are the 1939 grid cropped by 3 columns
# (same CRS and resolution, different extent).
# Run with: Rscript dev/issue-checks/issue-07.R
cat("Issue #7: ex_by_date()/bg_by_date(): first year's cell IDs are reused for all years, giving wrong background values when grids differ\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")
future::plan("sequential")

folders <- c("1939", "1940")
root <- file.path(tempdir(), "issue07_layers")
unlink(root, recursive = TRUE)
dir.create(root)
tmpl <- terra::rast(nrows = 30, ncols = 30, xmin = 0, xmax = 30, ymin = 0, ymax = 30,
                    crs = "EPSG:4326")
for (i in seq_along(folders)) {
  dir.create(file.path(root, folders[i]))
  r1 <- tmpl
  terra::values(r1) <- i * 10000 + seq_len(terra::ncell(tmpl))
  r2 <- tmpl
  terra::values(r2) <- i * 10000 + 5000 + seq_len(terra::ncell(tmpl))
  if (folders[i] == "1940") {
    r1 <- terra::crop(r1, terra::ext(3, 30, 0, 30))
    r2 <- terra::crop(r2, terra::ext(3, 30, 0, 30))
  }
  names(r1) <- "bio_01"
  names(r2) <- "bio_12"
  terra::writeRaster(r1, file.path(root, folders[i], "bio_01.tif"))
  terra::writeRaster(r2, file.path(root, folders[i], "bio_12.tif"))
}
occs <- data.frame(lon = c(10.5, 15.5, 20.5, 12.5), lat = c(10.5, 15.5, 20.5, 22.5),
                   date = c("1939", "1939", "1940", "1940"))
set.seed(7)
err <- NULL
bg <- tryCatch({
  # Since #25, sp_temporal_data() may already reject the differing grid
  st <- sp_temporal_data(occs = occs, longitude = "lon", latitude = "lat",
                         sp_date_var = "date", occ_date_format = "y",
                         layers_date_format = "y", layers_by_date_dir = root,
                         layers_ext = "*.tif$")
  ex <- suppressMessages(ex_by_date(st, train_prop = 0.7))
  suppressMessages(bg_by_date(ex, buffer_ngbs = 2, n_bg = 40))
}, error = function(e) {
  err <<- conditionMessage(e)
  NULL
})

if (is.null(bg)) {
  informative <- grepl("extent|resolution|grid|geometr", err, ignore.case = TRUE)
  got <- paste0("error: ", err)
  ok <- informative
} else {
  env_bg <- bg$env_bg
  folder <- basename(env_bg$ID_YEAR)
  at_xy <- rep(NA_real_, nrow(env_bg))
  for (f in unique(folder)) {
    rows <- which(folder == f)
    lyr <- terra::rast(file.path(root, f, "bio_01.tif"))
    at_xy[rows] <- terra::extract(lyr, as.matrix(env_bg[rows, c("lon", "lat")]))[, 1]
  }
  mism <- is.na(at_xy) | env_bg$bio_01 != at_xy
  got <- paste0("no error; background bio_01 differs from the value at its reported",
                " coordinates in ", sum(mism), " of ", nrow(env_bg), " rows (",
                paste(names(table(folder[mism])), table(folder[mism]), sep = ": ",
                      collapse = ", "), ")")
  ok <- !any(mism)
}

cat("EXPECTED: an informative error that the year layers do not share one grid,",
    "or background values that equal the values at their reported coordinates\n")
cat("GOT:     ", got, "\n")
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
