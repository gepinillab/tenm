# Issue #28: Columns are addressed by position (ex_by_date, tdf2swd,
# bg_by_date); env_data includes trian_test.
# Synthetic month folders "2000-01" and "2000-02" (layers_date_format = "ym"),
# two variables. The temporal_df columns are reordered (date first) and an
# extra column is added, as a user or an earlier tenm version might do.
# Run with: Rscript dev/issue-checks/issue-28.R
cat("Issue #28: Columns are addressed by position (ex_by_date, tdf2swd, bg_by_date); env_data includes trian_test\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")
future::plan("sequential")

root <- file.path(tempdir(), "issue28_layers")
unlink(root, recursive = TRUE)
dir.create(root)
tmpl <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20,
                    crs = "EPSG:4326")
for (i in 1:2) {
  f <- file.path(root, sprintf("2000-%02d", i))
  dir.create(f)
  r1 <- tmpl
  terra::values(r1) <- i * 1000 + seq_len(terra::ncell(tmpl))
  names(r1) <- "bio_01"
  r2 <- tmpl
  terra::values(r2) <- i * 1000 + 500 + seq_len(terra::ncell(tmpl))
  names(r2) <- "bio_12"
  terra::writeRaster(r1, file.path(f, "bio_01.tif"))
  terra::writeRaster(r2, file.path(f, "bio_12.tif"))
}
set.seed(28)
occ <- data.frame(lon = stats::runif(12, 1, 19), lat = stats::runif(12, 1, 19),
                  date = rep(c("2000-01", "2000-02"), each = 6))
st <- suppressMessages(sp_temporal_data(occ, "lon", "lat", "date", "ym", "ym", root))
tdf <- st$temporal_df
st$temporal_df <- data.frame(tdf[, c("date", "lon", "lat", "layer_dates", "layers_path")],
                             source = "survey")

run <- function(expr) tryCatch(suppressMessages(suppressWarnings(expr)),
                               error = function(e) conditionMessage(e))
ex <- run(ex_by_date(st, train_prop = 0.7))
ok_ex <- is.list(ex)
env_names <- if (ok_ex) names(ex$env_data) else ex
swd <- if (ok_ex) run(tdf2swd(ex, sp_name = "sp")) else "not run"
bg <- if (ok_ex) run(bg_by_date(ex, buffer_ngbs = 2, n_bg = 30)) else "not run"
bg_swd <- if (is.list(bg)) run(tdf2swd(bg)) else "not run"
# Each record's bio_01 must be the value of its own folder at its coordinates
env_values_ok <- ok_ex && is.data.frame(swd) && all(c("lon", "lat", "bio_01") %in% names(swd)) &&
  all(vapply(seq_len(nrow(swd)), function(k) {
    r <- terra::rast(file.path(ex$temporal_df$layers_path[k], "bio_01.tif"))
    swd$bio_01[k] == terra::extract(r, as.matrix(swd[k, c("lon", "lat")]))[1, 1]
  }, logical(1)))

nm <- function(x) if (is.data.frame(x)) paste(names(x), collapse = ",") else paste("error:", x)
ok <- ok_ex && identical(env_names, c("bio_01", "bio_12")) &&
  is.data.frame(swd) && identical(names(swd), c("sp_name", "lon", "lat", "date", "bio_01", "bio_12")) &&
  env_values_ok &&
  is.data.frame(bg_swd) && identical(names(bg_swd), c("sp_name", "lon", "lat", "year", "bio_01", "bio_12")) &&
  !anyNA(bg_swd$year) && all(bg_swd$year == 2000)

cat("EXPECTED: env_data has only bio_01, bio_12 (no trian_test); tdf2swd(ex) gives",
    "sp_name,lon,lat,date,bio_01,bio_12 with the values of each record's folder;",
    "tdf2swd(bg) gives sp_name,lon,lat,year,bio_01,bio_12 with year 2000 (from the layer date) and no NA\n")
cat("GOT:      env_data:", paste(env_names, collapse = ","), "| tdf2swd(ex):", nm(swd),
    "| values match =", env_values_ok, "| tdf2swd(bg):", nm(bg_swd),
    if (is.data.frame(bg_swd)) paste("| year:", paste(unique(bg_swd$year), collapse = ","),
                                     "| NA years:", sum(is.na(bg_swd$year))), "\n")
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
