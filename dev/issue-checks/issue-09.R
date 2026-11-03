# Issue #9: clean_dup(): by_mask = TRUE with a non-SpatRaster mask silently
# falls back to distance thinning (README affected).
# Uses the bundled abronia data and the 1939 bio_01 layer as the mask.
# Run with: Rscript dev/issue-checks/issue-09.R
cat("Issue #9: clean_dup(): by_mask = TRUE with a non-SpatRaster mask silently falls back to distance thinning (README affected)\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

data("abronia", package = "tenm")
d <- system.file("extdata/bio", package = "tenm")
tenm_mask <- terra::rast(file.path(d, "1939/bio_01.tif"))
cat("class(tenm_mask[1]):", class(tenm_mask[1]), "\n")

# The README call: raster_mask = tenm_mask[1] (a data.frame, not a SpatRaster)
err <- NULL
res <- tryCatch(
  clean_dup(abronia, "decimalLongitude", "decimalLatitude", by_mask = TRUE,
            threshold = terra::res(tenm_mask)[1], raster_mask = tenm_mask[1],
            progress = FALSE),
  error = function(e) {
    err <<- conditionMessage(e)
    NULL
  })
# by_mask = TRUE with no mask at all
err_null <- NULL
res_null <- tryCatch(
  clean_dup(abronia, "decimalLongitude", "decimalLatitude", by_mask = TRUE,
            threshold = 0.1, raster_mask = NULL, progress = FALSE),
  error = function(e) {
    err_null <<- conditionMessage(e)
    NULL
  })
by_mask_rows <- nrow(clean_dup(abronia, "decimalLongitude", "decimalLatitude",
                               by_mask = TRUE, raster_mask = tenm_mask,
                               progress = FALSE))

describe <- function(r, e) {
  if (is.null(e)) paste0("no error, ", nrow(r), " rows returned (distance thinning)")
  else paste0("error: ", e)
}
ok <- !is.null(err) && grepl("SpatRaster", err) &&
  !is.null(err_null) && grepl("SpatRaster", err_null)

cat("EXPECTED: an error saying raster_mask must be a SpatRaster, both for",
    "raster_mask = tenm_mask[1] and raster_mask = NULL (a real mask gives",
    by_mask_rows, "rows)\n")
cat("GOT:      tenm_mask[1] ->", describe(res, err), "\n")
cat("          NULL         ->", describe(res_null, err_null), "\n")
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
