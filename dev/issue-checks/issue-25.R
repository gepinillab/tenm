# Issue #25: No validation of CRS, extent, resolution or layer completeness
# across years.
# Grids that differ across years are already rejected by ex_by_date() (#7,
# see issue-07.R). This script checks the other cases with synthetic layers:
#   (a) one year folder lacks a variable;
#   (b) the layers have no CRS (a 0-1000 extent, so terra does not assume lon/lat);
#   (c) projected layers (EPSG:3857) with lon/lat occurrences.
# Run with: Rscript dev/issue-checks/issue-25.R
cat("Issue #25: No validation of CRS, extent, resolution or layer completeness across years\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")
future::plan("sequential")

make_layers <- function(name, crs, xmin, xmax, drop_in = NULL) {
  root <- file.path(tempdir(), name)
  unlink(root, recursive = TRUE)
  dir.create(root)
  for (f in c("2000", "2001")) {
    dir.create(file.path(root, f))
    r <- terra::rast(nrows = 10, ncols = 10, xmin = xmin, xmax = xmax, ymin = xmin,
                     ymax = xmax, crs = crs, vals = 1:100)
    names(r) <- "bio_01"
    terra::writeRaster(r, file.path(root, f, "bio_01.tif"))
    if (!identical(f, drop_in)) {
      names(r) <- "bio_12"
      terra::writeRaster(r, file.path(root, f, "bio_12.tif"))
    }
  }
  root
}
occ <- data.frame(lon = c(2.5, 5.5, 7.5, 3.5), lat = c(2.5, 5.5, 7.5, 6.5),
                  year = c(2000, 2000, 2001, 2001))

# Run sp_temporal_data() and ex_by_date(); collect warnings and the error
run <- function(root) {
  w <- character()
  err <- NULL
  ex <- withCallingHandlers(
    tryCatch({
      st <- sp_temporal_data(occ, "lon", "lat", "year", "y", "y", root)
      suppressMessages(ex_by_date(st, train_prop = 0.7))
    }, error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  list(ex = ex, w = paste(w, collapse = " | "), err = err)
}
describe <- function(r) {
  out <- if (!is.null(r$err)) paste("error:", r$err) else {
    nas <- colSums(is.na(r$ex$env_data))
    paste0("no error; NA values in env_data: ",
           paste(names(nas), nas, sep = " = ", collapse = ", "))
  }
  if (nzchar(r$w)) out <- paste0(out, "; warnings: ", r$w)
  gsub("\n", " ", out)
}

a <- run(make_layers("issue25_missing", "EPSG:4326", 0, 10, drop_in = "2001"))
b <- run(make_layers("issue25_nocrs", "", 0, 1000))
c3857 <- run(make_layers("issue25_projected", "EPSG:3857", 1e6, 2e6))

ok_a <- !is.null(a$err) && grepl("bio_12", a$err) && grepl("2001", a$err)
ok_b <- grepl("no CRS", b$w)
ok_c <- grepl("outside", c3857$w) && grepl("CRS", c3857$w)

cat("EXPECTED: (a) an error naming folder 2001 and the missing layer bio_12;",
    "(b) a warning that the layers have no CRS;",
    "(c) a warning that the occurrences fall outside the layers, suggesting a CRS check\n")
cat("GOT:      (a)", describe(a), "\n")
cat("          (b)", describe(b), "\n")
cat("          (c)", describe(c3857), "\n")
cat("RESULT:", if (ok_a && ok_b && ok_c) "PASS" else "FAIL", "\n")
