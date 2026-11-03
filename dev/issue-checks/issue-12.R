# Issue #12: bg_by_date(buffer_distance = ...): sample() on a length-1 vector
# picks a random wrong cell.
# Synthetic 1-degree layers written to tempdir(); value = cell number + 10000.
# Run with: Rscript dev/issue-checks/issue-12.R
cat("Issue #12: bg_by_date(buffer_distance = ...): sample() on a length-1 vector picks a random wrong cell\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")
future::plan("sequential")

root <- file.path(tempdir(), "issue12_layers")
unlink(root, recursive = TRUE)
dir.create(file.path(root, "1939"), recursive = TRUE)
tmpl <- terra::rast(nrows = 30, ncols = 30, xmin = 0, xmax = 30, ymin = 0, ymax = 30,
                    crs = "EPSG:4326")
r1 <- tmpl
terra::values(r1) <- 10000 + seq_len(terra::ncell(tmpl))
names(r1) <- "bio_01"
r2 <- tmpl
terra::values(r2) <- 15000 + seq_len(terra::ncell(tmpl))
names(r2) <- "bio_12"
terra::writeRaster(r1, file.path(root, "1939", "bio_01.tif"))
terra::writeRaster(r2, file.path(root, "1939", "bio_12.tif"))

run_bg <- function(occs, buffer_distance, n_bg, seed) {
  st <- sp_temporal_data(occs = occs, longitude = "lon", latitude = "lat",
                         sp_date_var = "date", occ_date_format = "y",
                         layers_date_format = "y", layers_by_date_dir = root,
                         layers_ext = "*.tif$")
  set.seed(seed)
  ex <- suppressMessages(ex_by_date(st, train_prop = 0.7))
  tryCatch(suppressMessages(bg_by_date(ex, buffer_ngbs = NULL,
                                       buffer_distance = buffer_distance,
                                       n_bg = n_bg))$env_bg,
           error = function(e) conditionMessage(e))
}

# A: one occurrence at a cell centre; a 1 km buffer covers only that cell.
occ_a <- data.frame(lon = 15.5, lat = 15.5, date = "1939")
cell_a <- terra::cellFromXY(tmpl, cbind(15.5, 15.5))
got_a <- vapply(1:5, function(seed) {
  bg <- run_bg(occ_a, buffer_distance = 1000, n_bg = 10, seed = seed)
  if (is.character(bg)) return(paste0("error: ", bg))
  paste0("cells ", paste(terra::cellFromXY(tmpl, as.matrix(bg[, c("lon", "lat")])),
                         collapse = ","), " (bio_01 ", paste(bg$bio_01, collapse = ","), ")")
}, character(1))
ok_a <- all(got_a == paste0("cells ", cell_a, " (bio_01 ", 10000 + cell_a, ")"))

# B: three occurrences in adjacent cells; 150 km buffers overlap, and n_bg is
# large enough to take every cell in the buffers.
occ_b <- data.frame(lon = c(10.5, 11.5, 12.5), lat = 10.5, date = "1939")
buf <- terra::buffer(terra::vect(as.matrix(occ_b[, 1:2]), crs = "EPSG:4326"), width = 150000)
n_unique_b <- length(unique(terra::cells(tmpl, buf)[, 2]))
bg_b <- run_bg(occ_b, buffer_distance = 150000, n_bg = 1000, seed = 1)
if (is.character(bg_b)) {
  got_b <- paste0("error: ", bg_b)
  ok_b <- FALSE
} else {
  dup_b <- sum(duplicated(bg_b[, c("lon", "lat")]))
  got_b <- paste0(nrow(bg_b), " rows, ", dup_b, " duplicated")
  ok_b <- nrow(bg_b) == n_unique_b && dup_b == 0
}

# C: cells2samp() when the neighbourhood holds a single cell (1 x 2 raster,
# occurrence in cell 1, buffer_ngbs = 1: the only neighbour is cell 2).
mask_c <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1,
                      crs = "EPSG:4326")
terra::values(mask_c) <- 1:2
got_c <- vapply(1:20, function(seed) {
  set.seed(seed)
  tryCatch(cells2samp(data = data.frame(lon = 0.5, lat = 0.5), longitude = "lon",
                      latitude = "lat", buffer_ngbs = 1, raster_mask = mask_c,
                      n_bg = 1, progress = FALSE),
           error = function(e) NA_real_)
}, numeric(1))
ok_c <- all(!is.na(got_c) & got_c == 2)

cat("EXPECTED: A) every seed returns only cell", cell_a, "(bio_01", 10000 + cell_a, ");",
    "B)", n_unique_b, "rows, 0 duplicated; C) cells2samp() returns cell 2 for all 20 seeds\n")
cat("GOT:      A)", paste(unique(got_a), collapse = " | "), "; B)", got_b,
    "; C) cells returned:", paste(names(table(got_c, useNA = "ifany")),
                                  table(got_c, useNA = "ifany"), sep = " x", collapse = ", "),
    if (anyNA(got_c)) "(NA = error)",
    "\n")
cat("RESULT:", if (ok_a && ok_b && ok_c) "PASS" else "FAIL", "\n")
