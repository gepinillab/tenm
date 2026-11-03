# Issue #11: clean_dup(): points outside the raster mask are silently
# collapsed or dropped.
# Synthetic 30 x 30 mask (0-30 lon/lat); three points far outside its
# extent and one inside.
# Run with: Rscript dev/issue-checks/issue-11.R
cat("Issue #11: clean_dup(): points outside the raster mask are silently collapsed or dropped\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

mask <- terra::rast(nrows = 30, ncols = 30, xmin = 0, xmax = 30, ymin = 0, ymax = 30,
                    crs = "EPSG:4326", vals = 1)
out <- data.frame(id = 1:4, lon = c(-120, -121, -122, 10.5), lat = c(40, 41, 42, 10.5))

run <- function(n_ngbs) {
  w <- character()
  res <- withCallingHandlers(
    clean_dup(out, "lon", "lat", by_mask = TRUE, raster_mask = mask,
              n_ngbs = n_ngbs, progress = FALSE),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  list(ids = res$id, warned = any(grepl("3 ", w) & grepl("outside", w)), w = w)
}
r0 <- run(0)
r1 <- run(1)
fmt <- function(r) paste0("keeps ids ", paste(r$ids, collapse = ", "), "; ",
                          if (r$warned) "warning names the 3 records outside the mask"
                          else if (length(r$w)) paste("warning:", r$w[1]) else "no warning")
ok <- identical(as.integer(r0$ids), 4L) && identical(as.integer(r1$ids), 4L) &&
  r0$warned && r1$warned

cat("EXPECTED: n_ngbs = 0 and n_ngbs = 1 both keep id 4 only and warn that 3",
    "records fall outside raster_mask\n")
cat("GOT:      n_ngbs = 0:", fmt(r0), "\n")
cat("          n_ngbs = 1:", fmt(r1), "\n")
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
