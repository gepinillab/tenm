# Issue #36: Global options and RNG state are changed without being restored
# (CRAN policy).
# pROC()'s RNG handling was fixed in #13 (see issue-13.R). This script checks
# what remains: the parallel branch of ellipsoid_selection() (more than 100
# models, parallel = TRUE) sets options(future.rng.onMisuse, future.globals.maxSize).
# Synthetic data: 15 variables taken by 2 = 105 models, 2 workers.
# It also checks that pROC() with rseed = TRUE leaves .Random.seed unchanged.
# Run with: Rscript dev/issue-checks/issue-36.R
cat("Issue #36: Global options and RNG state are changed without being restored (CRAN policy)\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

set.seed(36)
vars <- paste0("v", 1:15)
mk <- function(n) as.data.frame(stats::setNames(lapply(vars, function(v) stats::rnorm(n)), vars))
env_train <- mk(40)
env_test <- mk(15)
env_bg <- mk(300)

options(future.rng.onMisuse = NULL, future.globals.maxSize = NULL)
future::plan("sequential")
onames <- c("future.rng.onMisuse", "future.globals.maxSize")
before <- stats::setNames(lapply(onames, getOption), onames)
invisible(utils::capture.output(
  sel <- ellipsoid_selection(env_train, env_test, env_vars = vars, nvarstest = 2,
                             level = 0.95, mve = FALSE, env_bg = env_bg,
                             omr_criteria = 0.1, parallel = TRUE, ncores = 2,
                             proc = FALSE)))
after <- stats::setNames(lapply(onames, getOption), onames)
plan_ok <- inherits(future::plan(), "sequential")
opts_ok <- identical(before, after)
fmt <- function(o) paste(names(o), vapply(o, function(x) if (is.null(x)) "NULL" else
  format(x), character(1)), sep = " = ", collapse = ", ")

set.seed(1)
seed0 <- .Random.seed
invisible(pROC(stats::runif(500), stats::runif(30, 0.3, 1), n_iter = 20, rseed = TRUE))
rng_ok <- identical(seed0, .Random.seed)

cat("EXPECTED: after a parallel ellipsoid_selection() (", nrow(sel), "models), the options",
    "future.rng.onMisuse and future.globals.maxSize and the future plan are as before;",
    "pROC(rseed = TRUE) leaves .Random.seed unchanged\n")
cat("GOT:      options before: ", fmt(before), " | after: ", fmt(after),
    " | plan restored = ", plan_ok, " | .Random.seed unchanged by pROC() = ", rng_ok, "\n", sep = "")
cat("RESULT:", if (opts_ok && plan_ok && rng_ok) "PASS" else "FAIL", "\n")
