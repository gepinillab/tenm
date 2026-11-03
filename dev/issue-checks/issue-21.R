# Issue #21: clean_dup(): NA latitude and coordinate-only input cause errors.
# Uses the bundled abronia data and the 1939 bio_01 layer as the mask.
# Run with: Rscript dev/issue-checks/issue-21.R
cat("Issue #21: clean_dup(): NA latitude and coordinate-only input cause errors\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

data("abronia", package = "tenm")
d <- system.file("extdata/bio", package = "tenm")
tenm_mask <- terra::rast(file.path(d, "1939/bio_01.tif"))
lon <- "decimalLongitude"
lat <- "decimalLatitude"

try_it <- function(expr) tryCatch(suppressMessages(expr), error = function(e) {
  structure(list(msg = conditionMessage(e)), class = "failed")
})
describe <- function(r) {
  if (inherits(r, "failed")) paste("error:", gsub("\n", " ", r$msg))
  else paste0(nrow(r), " rows, columns: ", paste(names(r), collapse = ", "))
}

# (a) one NA latitude, distance thinning
x <- abronia[1:5, ]
x$decimalLatitude[2] <- NA
a <- try_it(clean_dup(x, lon, lat, threshold = 0.1))
# (b) coordinate-only input, by_mask
b <- try_it(clean_dup(abronia[, c(lon, lat)], lon, lat, by_mask = TRUE,
                      raster_mask = tenm_mask, progress = FALSE))
full <- clean_dup(abronia, lon, lat, by_mask = TRUE, raster_mask = tenm_mask,
                  progress = FALSE)
# (c) by_mask keeps the input column order (lon, lat, year, ...)
y <- abronia[, c(lon, lat, "year", "species")]
cc <- try_it(clean_dup(y, lon, lat, by_mask = TRUE, raster_mask = tenm_mask,
                       progress = FALSE))

ok_a <- !inherits(a, "failed") && !anyNA(a[, c(lon, lat)]) && nrow(a) >= 1
ok_b <- !inherits(b, "failed") && nrow(b) == nrow(full) &&
  isTRUE(all.equal(b[, c(lon, lat)], full[, c(lon, lat)], check.attributes = FALSE))
ok_c <- !inherits(cc, "failed") && identical(names(cc), names(y))

cat("EXPECTED: (a) the NA-latitude record is dropped, no error;",
    "(b) coordinate-only input gives the same", nrow(full), "records as the full data;",
    "(c) by_mask output keeps the input column order", paste(names(y), collapse = ", "), "\n")
cat("GOT:      (a)", describe(a), "\n")
cat("          (b)", describe(b), "\n")
cat("          (c)", describe(cc), "\n")
cat("RESULT:", if (ok_a && ok_b && ok_c) "PASS" else "FAIL", "\n")
