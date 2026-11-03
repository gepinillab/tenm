# Issue #6: ex_by_date(): train/test labels land on the wrong rows when folder
# order differs from date order.
# Synthetic layers in "my" folders written to tempdir().
# Run with: Rscript dev/issue-checks/issue-06.R
cat("Issue #6: ex_by_date(): train/test labels land on the wrong rows when folder order differs from date order\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")
future::plan("sequential")

# Folder (alphabetical) order: 01-1940, 03-1939, 12-1939.
# Date order:                  03-1939, 12-1939, 01-1940.
folders <- c("12-1939", "01-1940", "03-1939")
n_occ <- c(1, 2, 4)
train_prop <- 0.7
root <- file.path(tempdir(), "issue06_layers")
unlink(root, recursive = TRUE)
dir.create(root)
tmpl <- terra::rast(nrows = 30, ncols = 30, xmin = 0, xmax = 30, ymin = 0, ymax = 30,
                    crs = "EPSG:4326")
for (i in seq_along(folders)) {
  dir.create(file.path(root, folders[i]))
  r1 <- tmpl
  terra::values(r1) <- i * 10000 + seq_len(terra::ncell(tmpl))
  names(r1) <- "bio_01"
  r2 <- tmpl
  terra::values(r2) <- i * 10000 + 5000 + seq_len(terra::ncell(tmpl))
  names(r2) <- "bio_12"
  terra::writeRaster(r1, file.path(root, folders[i], "bio_01.tif"))
  terra::writeRaster(r2, file.path(root, folders[i], "bio_12.tif"))
}
occs <- data.frame(lon = c(5.5, 10.5, 12.5, 20.5, 22.5, 24.5, 26.5),
                   lat = c(5.5, 10.5, 12.5, 20.5, 22.5, 24.5, 26.5),
                   date = rep(folders, n_occ))
st <- sp_temporal_data(occs = occs, longitude = "lon", latitude = "lat",
                       sp_date_var = "date", occ_date_format = "my",
                       layers_date_format = "my", layers_by_date_dir = root,
                       layers_ext = "*.tif$")

# Expected labels per folder: 1 record -> 1 Train; 2 records -> 1 Train +
# 1 Test; n > 3 records -> ceiling(n * train_prop) Train, rest Test.
expected_train <- c(1, 1, ceiling(4 * train_prop))
results <- list()
for (seed in 1:5) {
  set.seed(seed)
  ex <- suppressMessages(ex_by_date(st, train_prop = train_prop))
  tdf <- ex$temporal_df
  folder <- basename(tdf$layers_path)
  # occurrence values must come from the record's own folder
  values_ok <- all(tdf$bio_01 == match(folder, folders) * 10000 +
                     terra::cellFromXY(tmpl, as.matrix(tdf[, c("lon", "lat")])))
  n_train <- vapply(folders, function(f) sum(tdf$trian_test[folder == f] == "Train"), numeric(1))
  n_test <- vapply(folders, function(f) sum(tdf$trian_test[folder == f] == "Test"), numeric(1))
  results[[seed]] <- data.frame(seed = seed,
                                labels = paste0(folders, ": ", n_train, " Train/", n_test, " Test",
                                                collapse = ", "),
                                ok = values_ok && all(n_train == expected_train) &&
                                  all(n_train + n_test == n_occ))
}
res <- do.call(rbind, results)

cat("EXPECTED: for every seed,",
    paste0(folders, ": ", expected_train, " Train/", n_occ - expected_train, " Test",
           collapse = ", "), "\n")
cat("GOT:\n")
print(res[, c("seed", "labels")], row.names = FALSE, right = FALSE)
cat("RESULT:", if (all(res$ok)) "PASS" else "FAIL", "\n")
