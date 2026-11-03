# Issue #23: sp_temporal_data(): occurrences with no matching layer date are
# silently matched to the nearest year.
# Synthetic year folders 1939, 1941 and 2016. This script checks that
# non-exact matches, unparseable dates, a wrong layers_date_format and a
# missing directory are reported. It does not check a matching policy
# (max_date_gap / on_unmatched); the nearest-date rule itself is unchanged.
# Run with: Rscript dev/issue-checks/issue-23.R
cat("Issue #23: sp_temporal_data(): occurrences with no matching layer date are silently matched to the nearest year\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

root <- file.path(tempdir(), "issue23_layers")
unlink(root, recursive = TRUE)
dir.create(root)
tmpl <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                    crs = "EPSG:4326", vals = 1)
names(tmpl) <- "bio_01"
for (f in c("1939", "1941", "2016")) {
  dir.create(file.path(root, f))
  terra::writeRaster(tmpl, file.path(root, f, "bio_01.tif"))
}
occ <- data.frame(lon = 5.5, lat = 5.5, year = c(1945, 1939, 1890, 2050, 2016, NA))

# Run a call and collect its messages, warnings and error
capture <- function(expr) {
  msgs <- character()
  err <- NULL
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(val = val, msgs = paste(trimws(msgs), collapse = " | "), err = err)
}

a <- capture(sp_temporal_data(occ, "lon", "lat", "year", "y", "y", root))
b <- capture(sp_temporal_data(occ, "lon", "lat", "year", "y", "dmy", root))
d <- capture(sp_temporal_data(occ, "lon", "lat", "year", "y", "y", file.path(root, "no_such_dir")))

ok_a <- !is.null(a$val) && grepl("2 of 5", a$msgs) && grepl("outside", a$msgs) &&
  grepl("1 record", a$msgs)
ok_b <- !is.null(b$err) && grepl("layers_date_format", b$err)
ok_d <- !is.null(d$err) && grepl("does not exist", d$err)

cat("EXPECTED: (a) 5 parseable dates, 2 exact: a report that 2 of 5 match exactly,",
    "that 1890 and 2050 are outside the layer dates, and that 1 record with a missing date was removed;",
    "(b) layers_date_format = 'dmy' -> an error naming layers_date_format;",
    "(c) a missing directory -> an error saying it does not exist\n")
cat("GOT:      (a)", if (is.null(a$val)) paste("error:", a$err) else
  paste0(nrow(a$val$temporal_df), " rows; messages: ", if (nzchar(a$msgs)) a$msgs else "none"), "\n")
cat("          (b)", if (is.null(b$err)) "no error" else paste("error:", b$err), "\n")
cat("          (c)", if (is.null(d$err)) "no error" else paste("error:", d$err), "\n")
cat("RESULT:", if (ok_a && ok_b && ok_d) "PASS" else "FAIL", "\n")
