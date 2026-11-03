# Issue #5: bg_by_date(): background is aligned to years by row position
# (wrong year, wrong per-year sample size).
# Synthetic 30 x 30 layers written to tempdir(). Each layer value encodes its
# year and cell (value = year_index * 10000 + cell [+ 5000 for bio_12]), so
# every background value can be traced back to the year layer it came from.
# Run with: Rscript dev/issue-checks/issue-05.R
cat("Issue #5: bg_by_date(): background is aligned to years by row position (wrong year, wrong per-year sample size)\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")
future::plan("sequential")

tmpl <- terra::rast(nrows = 30, ncols = 30, xmin = 0, xmax = 30, ymin = 0, ymax = 30,
                    crs = "EPSG:4326")
buffer_ngbs <- 2
n_bg <- 40
# 1 occurrence in the first year, 3 in the second; buffers do not overlap.
occ_xy <- data.frame(x = c(5.5, 24.5, 20.5, 24.5), y = c(5.5, 24.5, 24.5, 20.5),
                     year_index = c(1, 2, 2, 2))
expected_n <- c(ceiling(1 / 4 * n_bg), ceiling(3 / 4 * n_bg))   # 10 and 30

# Scenario A: "my" folders, where folder (alphabetical) order != date order.
# Scenario B: "y" folders, with the rows of temporal_df reversed.
scenarios <- list(
  A_my_folders = list(folders = c("12-1939", "01-1940"), fmt = "my", reverse = FALSE),
  B_rows_reversed = list(folders = c("1939", "1940"), fmt = "y", reverse = TRUE))

summary_rows <- list()
for (sc_name in names(scenarios)) {
  sc <- scenarios[[sc_name]]
  root <- file.path(tempdir(), paste0("issue05_", sc_name))
  unlink(root, recursive = TRUE)
  dir.create(root)
  for (i in seq_along(sc$folders)) {
    dir.create(file.path(root, sc$folders[i]))
    r1 <- tmpl
    terra::values(r1) <- i * 10000 + seq_len(terra::ncell(tmpl))
    names(r1) <- "bio_01"
    r2 <- tmpl
    terra::values(r2) <- i * 10000 + 5000 + seq_len(terra::ncell(tmpl))
    names(r2) <- "bio_12"
    terra::writeRaster(r1, file.path(root, sc$folders[i], "bio_01.tif"))
    terra::writeRaster(r2, file.path(root, sc$folders[i], "bio_12.tif"))
  }
  occs <- data.frame(lon = occ_xy$x, lat = occ_xy$y,
                     date = sc$folders[occ_xy$year_index])
  set.seed(5)
  st <- sp_temporal_data(occs = occs, longitude = "lon", latitude = "lat",
                         sp_date_var = "date", occ_date_format = sc$fmt,
                         layers_date_format = sc$fmt, layers_by_date_dir = root,
                         layers_ext = "*.tif$")
  ex <- suppressMessages(ex_by_date(st, train_prop = 0.7))
  if (sc$reverse) ex$temporal_df <- ex$temporal_df[rev(seq_len(nrow(ex$temporal_df))), ]
  bg <- suppressMessages(bg_by_date(ex, buffer_ngbs = buffer_ngbs, n_bg = n_bg))
  env_bg <- bg$env_bg

  folder <- basename(env_bg$ID_YEAR)
  year_index <- match(folder, sc$folders)
  cell <- terra::cellFromXY(tmpl, as.matrix(env_bg[, c("lon", "lat")]))
  rc <- terra::rowColFromCell(tmpl, cell)
  occ_rc <- terra::rowColFromCell(tmpl, terra::cellFromXY(tmpl, as.matrix(occ_xy[, c("x", "y")])))
  # background point is within buffer_ngbs cells of an occurrence of its own year
  in_own_buffer <- vapply(seq_len(nrow(env_bg)), function(j) {
    own <- occ_rc[occ_xy$year_index == year_index[j], , drop = FALSE]
    any(pmax(abs(own[, 1] - rc[j, 1]), abs(own[, 2] - rc[j, 2])) <= buffer_ngbs)
  }, logical(1))
  # values come from that year's layers, at the reported coordinates
  value_ok <- env_bg$bio_01 == year_index * 10000 + cell &
    env_bg$bio_12 == year_index * 10000 + 5000 + cell
  counts <- table(factor(folder, levels = sc$folders))
  summary_rows[[sc_name]] <- data.frame(
    scenario = sc_name,
    counts = paste(names(counts), counts, sep = "=", collapse = ", "),
    expected = paste(sc$folders, expected_n, sep = "=", collapse = ", "),
    in_own_buffer = sprintf("%d/%d", sum(in_own_buffer), nrow(env_bg)),
    value_matches_year = sprintf("%d/%d", sum(value_ok), nrow(env_bg)),
    ok = all(as.vector(counts) == expected_n) && all(in_own_buffer) && all(value_ok))
}
res <- do.call(rbind, summary_rows)
rownames(res) <- NULL

cat("EXPECTED: for both scenarios, per-year counts",
    paste(expected_n, collapse = " and "),
    "(proportional to occurrences), every background point inside its own year's buffer,",
    "and every value read from its own year's layers\n")
cat("GOT:\n")
print(res[, c("scenario", "counts", "in_own_buffer", "value_matches_year")], row.names = FALSE)
cat("RESULT:", if (all(res$ok)) "PASS" else "FAIL", "\n")
