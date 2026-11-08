# Issue #56: Models with no test records keep pval_bin = 1 (a binomial test
# on 0 trials). Synthetic data: the test table has records, but none of them
# is complete for one of the variables, so the model on that variable has no
# test records.
# Run with: Rscript dev/issue-checks/issue-56.R
cat("Issue #56: Models with no test records keep pval_bin = 1\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

set.seed(56)
env_train <- data.frame(bio_01 = rnorm(60, 20, 2), bio_12 = rnorm(60, 1000, 100))
env_train$bio_30 <- rnorm(60, 5, 1)
env_bg <- data.frame(bio_01 = runif(1000, 5, 45), bio_12 = runif(1000, 200, 2500),
                     bio_30 = runif(1000, 0, 10))
env_test <- data.frame(bio_01 = rnorm(20, 21, 3), bio_12 = rnorm(20, 1050, 150),
                       bio_30 = rnorm(20, 5, 1.5))
env_test$bio_12 <- NA_real_   # no test record has bio_12

cat("EXPECTED: omr_test and pval_bin are NA for the bio_12 model (no test",
    "records),\n          and pval_bin is a number for the model without bio_12\n")
# (a) ellipsoid_omr() called directly with an empty test table
direct <- ellipsoid_omr(env_data = env_train[, "bio_12", drop = FALSE],
                        env_test = env_test[0, "bio_12", drop = FALSE],
                        env_bg = env_bg[, "bio_12", drop = FALSE],
                        cf_level = 0.975, mve = FALSE, proc = FALSE)
cat("GOT (ellipsoid_omr, 0 test rows): om_rate_test =", direct$om_rate_test,
    "pval_bin =", direct$pval_bin, "\n")
# (b) through ellipsoid_selection()
sel <- NULL
invisible(utils::capture.output(
  sel <- ellipsoid_selection(env_train = env_train, env_test = env_test,
                             env_vars = c("bio_01", "bio_12", "bio_30"),
                             nvarstest = 2,
                             level = 0.975, mve = FALSE, env_bg = env_bg,
                             omr_criteria = 0.5, parallel = FALSE,
                             proc = FALSE)))
print(sel[, c("fitted_vars", "om_rate_test", "pval_bin", "rank_by_omr_train_test")])
no_test <- grepl("bio_12", sel$fitted_vars)
with_test <- !no_test
ok_a <- is.na(direct$pval_bin)
ok_b <- all(is.na(sel$pval_bin[no_test])) && all(!is.na(sel$pval_bin[with_test])) &&
  max(sel$rank_by_omr_train_test[with_test]) < min(sel$rank_by_omr_train_test[no_test])
cat("\nRESULT:", if (ok_a && ok_b) "PASS" else "FAIL", "\n")
