# Issue #4: ellipsoid_omr(): binomial-test p-value is off by one
# (P(X > k) instead of P(X >= k)).
# Synthetic data only. Run with: Rscript dev/issue-checks/issue-04.R
cat("Issue #4: ellipsoid_omr(): binomial-test p-value is off by one (P(X > k) instead of P(X >= k))\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

set.seed(4)
env_train <- data.frame(bio_01 = rnorm(60, 20, 3), bio_12 = rnorm(60, 1000, 150))
env_bg <- data.frame(bio_01 = runif(2000, 5, 40), bio_12 = runif(2000, 300, 2000))

# Case A: all 5 test points inside the ellipsoid (the case the issue reports).
env_test_a <- data.frame(bio_01 = c(19, 20, 21, 20, 20.5),
                         bio_12 = c(1000, 980, 1020, 1050, 990))
# Case B: 12 test points, some outside the ellipsoid.
env_test_b <- data.frame(bio_01 = c(rnorm(9, 20, 2), 35, 6, 38),
                         bio_12 = c(rnorm(9, 1000, 100), 1900, 350, 400))

check_case <- function(env_test) {
  r <- ellipsoid_omr(env_data = env_train, env_test = env_test, env_bg = env_bg,
                     cf_level = 0.975, mve = FALSE, proc = FALSE)
  n <- nrow(env_test)
  k <- round(n * (1 - r$om_rate_test))   # test points inside the ellipsoid
  ref <- stats::binom.test(k, n, p = r$bg_prevalence, alternative = "greater")$p.value
  c(k = k, n = n, prev = r$bg_prevalence, pval_bin = r$pval_bin, binom_test = ref)
}
a <- check_case(env_test_a)
b <- check_case(env_test_b)

ok <- isTRUE(all.equal(a[["pval_bin"]], a[["binom_test"]])) &&
  isTRUE(all.equal(b[["pval_bin"]], b[["binom_test"]]))

cat("EXPECTED: pval_bin equals binom.test(k, n, bg_prevalence, alternative = 'greater')$p.value:",
    sprintf("A (k=%d, n=%d, prev=%.4f) %.6f; B (k=%d, n=%d, prev=%.4f) %.6f",
            a[["k"]], a[["n"]], a[["prev"]], a[["binom_test"]],
            b[["k"]], b[["n"]], b[["prev"]], b[["binom_test"]]), "\n")
cat("GOT:      pval_bin:", sprintf("A %.6f; B %.6f", a[["pval_bin"]], b[["pval_bin"]]), "\n")
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
