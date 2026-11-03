# Issue #26: ellipsoid_selection(): missing test omission is set to 0, failed
# models vanish, and the output schema varies.
# Synthetic data only. Run with: Rscript dev/issue-checks/issue-26.R
cat("Issue #26: ellipsoid_selection(): missing test omission is set to 0, failed models vanish, and the output schema varies\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

set.seed(26)
env_train <- data.frame(bio_01 = rnorm(60, 20, 3), bio_12 = rnorm(60, 1000, 150),
                        bio_03 = rnorm(60, 50, 5))
# bio_03 has no test records; bio_99 = 2 * bio_01, so the (bio_01, bio_99)
# covariance matrix is singular and that model cannot be fitted.
env_test <- data.frame(bio_01 = rnorm(20, 20, 3), bio_12 = rnorm(20, 1000, 150),
                       bio_03 = NA_real_)
env_bg <- data.frame(bio_01 = runif(2000, 5, 40), bio_12 = runif(2000, 200, 2500),
                     bio_03 = runif(2000, 20, 80))
env_train$bio_99 <- 2 * env_train$bio_01
env_test$bio_99 <- 2 * env_test$bio_01
env_bg$bio_99 <- 2 * env_bg$bio_01
vars <- c("bio_01", "bio_12", "bio_03", "bio_99")

run <- function(omr_criteria) {
  out <- NULL
  invisible(utils::capture.output(
    out <- suppressWarnings(ellipsoid_selection(env_train, env_test, env_vars = vars,
                                                nvarstest = 2, level = 0.975, mve = FALSE,
                                                env_bg = env_bg, omr_criteria = omr_criteria,
                                                parallel = FALSE, proc = TRUE, proc_iter = 5)),
    type = "output"))
  out
}
sink(tempfile()); sink(type = "message", file(tempfile(), "w"))  # pROC() try() messages
sel <- run(0.1)    # some models pass omr_criteria
sel_none <- run(-1)  # no model passes omr_criteria
sink(type = "message"); sink()

no_test <- grepl("bio_03", sel$fitted_vars)
with_test <- !no_test & sel$fitted_vars != "bio_01,bio_99"
# A) models without test records: om_rate_test is NA and they rank after models with test data
rk <- sel$rank_by_omr_train_test
ok_a <- all(is.na(sel$om_rate_test[no_test])) &&
  max(which(with_test)) < min(which(no_test)) &&
  max(rk[with_test]) < min(rk[no_test])
# B) the singular model is reported, not dropped
failed_row <- sel[sel$fitted_vars == "bio_01,bio_99", , drop = FALSE]
ok_b <- nrow(failed_row) == 1 && "status" %in% names(sel) &&
  identical(failed_row$status, "fit_failed") && !is.na(failed_row$fit_error)
# C) the same columns whether or not a model passes omr_criteria
ok_c <- identical(names(sel), names(sel_none)) && "rank_omr_aucratio" %in% names(sel_none)

cat("EXPECTED: A) om_rate_test = NA for the 3 models with bio_03 (no test records), ranked after",
    "models with test data (rows and rank_by_omr_train_test); B) a 'fit_failed' row with its error for bio_01,bio_99;",
    "C) identical columns, including rank_omr_aucratio, when no model passes omr_criteria\n")
cat("GOT:      A) om_rate_test for bio_03 models =",
    paste(format(sel$om_rate_test[no_test]), collapse = ", "),
    "| rows of bio_03 models =", paste(which(no_test), collapse = ", "),
    "of", nrow(sel), "| their rank_by_omr_train_test =", paste(rk[no_test], collapse = ", "),
    "vs", paste(rk[with_test], collapse = ", "), "for models with test data\n")
cat("          B) bio_01,bio_99:", if (nrow(failed_row) == 0) "row missing" else
  paste0("status = ", failed_row$status, ", fit_error = \"", failed_row$fit_error, "\""), "\n")
cat("          C) columns only when some pass:",
    paste(setdiff(names(sel), names(sel_none)), collapse = ", "),
    "| columns only when none pass:", paste(setdiff(names(sel_none), names(sel)), collapse = ", "), "\n")
cat("RESULT:", if (ok_a && ok_b && ok_c) "PASS" else "FAIL", "\n")
