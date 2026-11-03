# Issue #10: clean_dup(): distance thinning over-thins chains of nearby records.
# Synthetic points on a line: a chain A-B-C (A-B and B-C closer than the
# threshold, A-C farther), plus a second chain of five points 0.6 apart.
# Run with: Rscript dev/issue-checks/issue-10.R
cat("Issue #10: clean_dup(): distance thinning over-thins chains of nearby records\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

threshold <- 1
chain3 <- data.frame(id = 1:3, x = c(0, 0.6, 1.2), y = 0)
chain5 <- data.frame(id = 1:5, x = c(0, 0.6, 1.2, 1.8, 2.4) + 10, y = 0)

# A thinned set is valid when no two kept records are within the threshold,
# and maximal when every dropped record is within the threshold of a kept one.
check_set <- function(dat, kept_ids) {
  d <- as.matrix(stats::dist(dat[, c("x", "y")]))
  k <- dat$id %in% kept_ids
  dk <- d[k, k, drop = FALSE]
  valid <- all(dk[upper.tri(dk)] >= threshold)
  maximal <- all(vapply(which(!k), function(i) any(d[i, k] < threshold), logical(1)))
  valid && maximal
}
k3 <- clean_dup(chain3, "x", "y", threshold = threshold)$id
k5 <- clean_dup(chain5, "x", "y", threshold = threshold)$id
ok <- identical(as.integer(k3), c(1L, 3L)) && check_set(chain3, k3) && check_set(chain5, k5)

cat("EXPECTED: chain of 3 keeps ids 1, 3; both results are valid and maximal",
    "(no kept pair closer than the threshold, every dropped record near a kept one)\n")
cat("GOT:      chain of 3 keeps ids", paste(k3, collapse = ", "),
    if (check_set(chain3, k3)) "(valid, maximal)" else "(not maximal)", "\n")
cat("          chain of 5 keeps ids", paste(k5, collapse = ", "),
    if (check_set(chain5, k5)) "(valid, maximal)" else "(not maximal)", "\n")
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
