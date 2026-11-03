# Issue #22: pROC(): proc_precision() errors when test data have fewer than
# 2 unique values.
# Synthetic suitability values only. Run with: Rscript dev/issue-checks/issue-22.R
cat("Issue #22: pROC(): proc_precision() errors when test data have fewer than 2 unique values\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

set.seed(1)
bgv <- runif(5000)
tv <- runif(30, .3, 1)
try_proc <- function(...) tryCatch(pROC(...)$pROC_summary, error = function(e) e)
show <- function(x) if (inherits(x, "error")) paste0("error \"", conditionMessage(x), "\"") else
  paste(names(x), format(x, digits = 4), sep = " = ", collapse = ", ")
finite <- function(x) !inherits(x, "error") && all(is.finite(x))

# A) all test values identical (1 unique value)
a <- try_proc(bgv, rep(0.5, 10), n_iter = 5, rseed = TRUE)
# B) values whose midpoint prints without a decimal point (0-100 scale):
#    mean(min background = 0, second test value = 60) = 30
b <- try_proc(0:100, c(40, 60, 80, 90), n_iter = 5, rseed = TRUE, sub_sample = FALSE)
# C) the usual case keeps working and gives the same result as before
#    (reference values from tenm 0.5.1.9001)
c_ <- try_proc(bgv, tv, n_iter = 20, rseed = TRUE)
c_ref <- c(Mean_AUC = 0.618458333333335, `Mean_pAUC_ratio_at_5%` = 1.22192439122705, P_value = 0)
# D) inside ellipsoid_omr(): identical test points give NA pROC columns
set.seed(22)
env_train <- data.frame(bio_01 = rnorm(60, 20, 3), bio_12 = rnorm(60, 1000, 150))
env_test <- data.frame(bio_01 = rep(20, 5), bio_12 = rep(1000, 5))
env_bg <- data.frame(bio_01 = runif(2000, 5, 40), bio_12 = runif(2000, 200, 2500))
d <- suppressWarnings(utils::capture.output(
  r <- ellipsoid_omr(env_train, env_test, env_bg, cf_level = 0.975, mve = FALSE,
                     proc = TRUE, proc_iter = 5), type = "message"))

ok <- finite(a) && finite(b) && finite(c_) &&
  isTRUE(all.equal(unname(c_), unname(c_ref))) &&
  !anyNA(r[, c("pval_proc", "env_bg_paucratio", "env_bg_auc")])

cat("EXPECTED: finite pROC summaries for A) 10 identical test values, B) integer 0-100",
    "suitabilities, C) the usual case (unchanged:", show(c_ref), ");",
    "D) ellipsoid_omr() with 5 identical test points has non-NA pROC columns\n")
cat("GOT:      A)", show(a), "\n")
cat("          B)", show(b), "\n")
cat("          C)", show(c_), "\n")
cat("          D) pval_proc =", r$pval_proc, "| env_bg_paucratio =", format(r$env_bg_paucratio, digits = 4),
    "| env_bg_auc =", format(r$env_bg_auc, digits = 4), "\n")
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
