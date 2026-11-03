# Issue #2: ellipsoid_selection(): proc_iter is passed by position and binds
# to mve, so mve and NoOfIteration are ignored.
# Synthetic data only. Run with: Rscript dev/issue-checks/issue-02.R
cat("Issue #2: ellipsoid_selection(): proc_iter is passed by position and binds to mve, so mve and NoOfIteration are ignored\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

# Training data with a cluster of outliers, so the minimum volume ellipsoid
# (mve = TRUE) and the classical covariance ellipsoid (mve = FALSE) differ.
set.seed(42)
env_train <- data.frame(bio_01 = c(rnorm(50, 20, 2), rnorm(10, 30, 6)),
                        bio_12 = c(rnorm(50, 1000, 100), rnorm(10, 1500, 300)))
env_test <- data.frame(bio_01 = rnorm(20, 21, 3),
                       bio_12 = rnorm(20, 1050, 150))
env_bg <- data.frame(bio_01 = runif(3000, 5, 45),
                     bio_12 = runif(3000, 200, 2500))

# Reference: ellipsoid_omr() called directly with named arguments.
direct_cov <- ellipsoid_omr(env_data = env_train, env_test = env_test,
                            env_bg = env_bg, cf_level = 0.975, mve = FALSE,
                            proc = FALSE)
set.seed(1)
direct_mve <- ellipsoid_omr(env_data = env_train, env_test = env_test,
                            env_bg = env_bg, cf_level = 0.975, mve = TRUE,
                            proc = FALSE)
cat("direct ellipsoid_omr(mve = FALSE) bg_prevalence:", direct_cov$bg_prevalence, "\n")
cat("direct ellipsoid_omr(mve = TRUE)  bg_prevalence:", direct_mve$bg_prevalence, "\n")

# Record which mve / proc_iter values ellipsoid_omr() receives, and which
# n_iter pROC() receives, when called from ellipsoid_selection().
rec <- new.env()
rec$omr <- list()
rec$proc_n_iter <- c()
invisible(suppressMessages(trace("ellipsoid_omr",
      tracer = quote(assign("omr", c(get("omr", envir = rec),
                                     list(list(mve = mve, proc_iter = proc_iter))),
                            envir = rec)),
      where = asNamespace("tenm"), print = FALSE)))
invisible(suppressMessages(trace("pROC",
      tracer = quote(assign("proc_n_iter", c(get("proc_n_iter", envir = rec), n_iter),
                            envir = rec)),
      where = asNamespace("tenm"), print = FALSE)))

set.seed(1)
invisible(utils::capture.output(
  sel <- ellipsoid_selection(env_train = env_train, env_test = env_test,
                             env_vars = c("bio_01", "bio_12"), nvarstest = 2,
                             level = 0.975, mve = FALSE, env_bg = env_bg,
                             omr_criteria = 0.1, parallel = FALSE,
                             proc = TRUE, proc_iter = 7, rseed = TRUE)))
suppressMessages(untrace("ellipsoid_omr", where = asNamespace("tenm")))
suppressMessages(untrace("pROC", where = asNamespace("tenm")))

got_mve <- rec$omr[[1]][["mve"]]
got_proc_iter <- rec$omr[[1]][["proc_iter"]]
got_n_iter <- rec$proc_n_iter[1]

ok_discriminates <- !isTRUE(all.equal(direct_cov$bg_prevalence, direct_mve$bg_prevalence))
ok_prev <- isTRUE(all.equal(sel$bg_prevalence[1], direct_cov$bg_prevalence))
ok_args <- isFALSE(got_mve) && isTRUE(all.equal(got_proc_iter, 7)) &&
  isTRUE(all.equal(got_n_iter, 7))

cat("\nEXPECTED: ellipsoid_selection(mve = FALSE, proc_iter = 7) passes mve = FALSE and proc_iter = 7",
    "to ellipsoid_omr(), pROC() runs n_iter = 7, and bg_prevalence equals the direct",
    "mve = FALSE value", format(direct_cov$bg_prevalence, digits = 6), "\n")
cat("GOT:      ellipsoid_omr() received mve =", got_mve, "| proc_iter =", got_proc_iter,
    "| pROC() n_iter =", got_n_iter,
    "| bg_prevalence =", format(sel$bg_prevalence[1], digits = 6), "\n")
if (!ok_discriminates) cat("NOTE: data do not separate mve = TRUE from mve = FALSE\n")
cat("RESULT:", if (ok_discriminates && ok_prev && ok_args) "PASS" else "FAIL", "\n")
