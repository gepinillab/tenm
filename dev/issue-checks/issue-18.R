# Issue #18: ex_by_date(): exact duplicate records produce list-columns that
# break later steps.
# Synthetic year folder 2014; the first record appears twice.
# Run with: Rscript dev/issue-checks/issue-18.R
cat("Issue #18: ex_by_date(): exact duplicate records produce list-columns that break later steps\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")
future::plan("sequential")

root <- file.path(tempdir(), "issue18_layers")
unlink(root, recursive = TRUE)
dir.create(root)
dir.create(file.path(root, "2014"))
r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                 crs = "EPSG:4326", vals = 1:100)
names(r) <- "bio_01"
terra::writeRaster(r, file.path(root, "2014", "bio_01.tif"))
r12 <- r * 10
names(r12) <- "bio_12"
terra::writeRaster(r12, file.path(root, "2014", "bio_12.tif"))

occ <- data.frame(lon = c(1.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
                  lat = c(1.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5), year = 2014)
st <- sp_temporal_data(occ, "lon", "lat", "year", "y", "y", root)
w <- character()
set.seed(18)
ex <- withCallingHandlers(suppressMessages(ex_by_date(st, train_prop = 0.7)),
                          warning = function(cnd) {
                            w <<- c(w, conditionMessage(cnd))
                            invokeRestart("muffleWarning")
                          })
cls <- vapply(ex$temporal_df[, c("bio_01", "bio_12")], function(x) class(x)[1], "")
expected_bio01 <- terra::extract(r, as.matrix(occ[, c("lon", "lat")]))[, 1]
ok <- all(cls == "numeric") && nrow(ex$temporal_df) == nrow(occ) &&
  isTRUE(all.equal(sort(ex$temporal_df$bio_01), sort(expected_bio01)))

cat("EXPECTED: bio_01 and bio_12 are numeric columns with", nrow(occ),
    "rows (one per record, duplicates included) and no warning\n")
cat("GOT:      classes", paste(names(cls), cls, sep = " = ", collapse = ", "), "; ",
    nrow(ex$temporal_df), "rows; warnings:", if (length(w)) paste(w, collapse = " | ") else "none", "\n")
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
