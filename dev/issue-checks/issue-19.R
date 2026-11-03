# Issue #19: cov_center(): the default vars = NULL selects zero columns and
# gives a misleading error.
# Synthetic data: three numeric variables, 60 rows.
# Run with: Rscript dev/issue-checks/issue-19.R
cat("Issue #19: cov_center(): the default vars = NULL selects zero columns and gives a misleading error\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

set.seed(19)
d <- data.frame(bio_01 = rnorm(60, 20, 2), bio_05 = rnorm(60, 30, 3), bio_12 = rnorm(60, 1000, 100))
run <- function(expr) tryCatch(expr, error = function(e) paste("error:", conditionMessage(e)))

set.seed(1)
default <- run(cov_center(d, mve = TRUE, level = 0.975))
set.seed(1)
named <- cov_center(d, mve = TRUE, level = 0.975, vars = names(d))
# Collinear variables (bio_05 = 2 * bio_01): singular covariance
sing <- run(cov_center(data.frame(a = d$bio_01, b = 2 * d$bio_01), mve = FALSE, level = 0.975))

ok_default <- is.list(default) && isTRUE(all.equal(default, named))
ok_sing <- is.character(sing) && grepl("singular", sing, ignore.case = TRUE)
cat("EXPECTED: vars = NULL uses all columns (same result as vars = names(data));",
    "collinear variables give an error that says the covariance matrix is singular\n")
cat("GOT:      vars = NULL:", if (is.list(default)) paste0("centroid of length ", length(default$centroid),
                                                         if (ok_default) " (same as vars = names(data))" else " (differs)")
                          else gsub("\n", " ", default), "\n")
cat("          collinear: ", if (is.character(sing)) gsub("\n", " ", sing) else "no error", "\n")
cat("RESULT:", if (ok_default && ok_sing) "PASS" else "FAIL", "\n")
