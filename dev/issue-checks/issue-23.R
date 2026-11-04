# Issue #23: sp_temporal_data(): occurrences with no matching layer date are
# silently matched to the nearest year.
# Synthetic year folders 1939, 1941 and 2016. Policy: each occurrence date is
# floored to the layer unit (its year for "y" folders) and matched exactly;
# on_unmatched = "drop" (default) removes records with no folder for their
# year and reports them, "nearest" keeps the old nearest-date matching and
# reports where records were moved, "error" stops. Also checks that
# clean_dup_by_date() then only compares records of the same year, and the
# input checks from the previous batch (bad layers_date_format, missing
# directory).
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
# All records in the same pixel. 1940-11-20 is nearer to the 1941 folder
# than to 1939, but its year has no folder.
occ <- data.frame(lon = 5.5, lat = 5.5,
                  date = c("1939-05-01", "1940-11-20", "1941-03-01", "1941-09-01",
                           "1890-01-01", "2050-01-01", "2016-06-30"))

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
yr <- function(x) format(as.Date(x), "%Y")

# (a) default: exact match on the year, unmatched records dropped and reported
a <- capture(sp_temporal_data(occ, "lon", "lat", "date", "ymd", "y", root))
a_df <- if (!is.null(a$val)) a$val$temporal_df else NULL
ok_a <- !is.null(a_df) && nrow(a_df) == 4 &&
  all(yr(a_df$date) == yr(a_df$layer_dates)) &&
  grepl("3 of 7", a$msgs) && grepl("1940", a$msgs) && grepl("1890", a$msgs) &&
  grepl("2050", a$msgs) && grepl("dropped", a$msgs)

# (b) clean_dup_by_date(): each layer folder only holds records of its year,
# so the 1939 and 1941 records in the same pixel are both kept, and the two
# 1941 records are reduced to one.
b <- if (!is.null(a$val)) capture(clean_dup_by_date(a$val, by_mask = TRUE, raster_mask = tmpl,
                                                    threshold = 1, n_ngbs = 0)) else list(val = NULL)
b_df <- if (!is.null(b$val)) b$val$temporal_df else NULL
ok_b <- !is.null(b_df) && nrow(b_df) == 3 &&
  setequal(yr(b_df$layer_dates), c("1939", "1941", "2016")) &&
  all(yr(b_df$date) == yr(b_df$layer_dates))

# (c) on_unmatched = "nearest": old matching, with the moves reported
c_ <- tryCatch(capture(sp_temporal_data(occ, "lon", "lat", "date", "ymd", "y", root,
                                        on_unmatched = "nearest")),
               error = function(e) list(val = NULL, msgs = "", err = conditionMessage(e)))
c_df <- if (!is.null(c_$val)) c_$val$temporal_df else NULL
ok_c <- !is.null(c_df) && nrow(c_df) == 7 &&
  yr(c_df$layer_dates[c_df$date == "1940-11-20"]) == "1941" &&
  grepl("1940 -> 1941", c_$msgs) && grepl("1890 -> 1939", c_$msgs)

# (d) on_unmatched = "error"
d <- tryCatch(capture(sp_temporal_data(occ, "lon", "lat", "date", "ymd", "y", root,
                                       on_unmatched = "error")),
              error = function(e) list(val = NULL, msgs = "", err = conditionMessage(e)))
ok_d <- !is.null(d$err) && grepl("no layer folder", d$err)

# (e) input checks
e1 <- capture(sp_temporal_data(occ, "lon", "lat", "date", "ymd", "dmy", root))
e2 <- capture(sp_temporal_data(occ, "lon", "lat", "date", "ymd", "y", file.path(root, "no_such_dir")))
ok_e <- !is.null(e1$err) && grepl("layers_date_format", e1$err) &&
  !is.null(e2$err) && grepl("does not exist", e2$err)

cat("EXPECTED: (a) default: 4 of 7 records kept (1939, 1941, 1941, 2016), each matched to the",
    "folder of its own year; a message that 3 of 7 were dropped, naming 1940, 1890 and 2050;",
    "(b) clean_dup_by_date() keeps 3 records (one per year 1939, 1941, 2016) and never mixes years;",
    "(c) on_unmatched = 'nearest': 7 records, 1940-11-20 -> 1941, moves reported ('1940 -> 1941',",
    "'1890 -> 1939'); (d) on_unmatched = 'error' stops; (e) bad layers_date_format and missing",
    "directory give explicit errors\n")
cat("GOT:      (a)", if (is.null(a_df)) paste("error:", a$err) else
  paste0(nrow(a_df), " rows, years ", paste(yr(a_df$date), collapse = ","), " -> layers ",
         paste(yr(a_df$layer_dates), collapse = ","), "; messages: ",
         if (nzchar(a$msgs)) a$msgs else "none"), "\n")
cat("          (b)", if (is.null(b_df)) "no result" else
  paste0(nrow(b_df), " rows, years ", paste(yr(b_df$date), collapse = ","), " -> layers ",
         paste(yr(b_df$layer_dates), collapse = ",")), "\n")
cat("          (c)", if (is.null(c_df)) paste("error:", c_$err) else
  paste0(nrow(c_df), " rows; messages: ", if (nzchar(c_$msgs)) c_$msgs else "none"), "\n")
cat("          (d)", if (is.null(d$err)) "no error" else paste("error:", d$err), "\n")
cat("          (e)", if (is.null(e1$err)) "no error" else e1$err, "|",
    if (is.null(e2$err)) "no error" else e2$err, "\n")
cat("RESULT:", if (ok_a && ok_b && ok_c && ok_d && ok_e) "PASS" else "FAIL", "\n")
