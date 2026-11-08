# Issue #54: ellipsoid_selection() fails with nvarstest = 1 because
# env_train[, var_comb] drops to a vector. Synthetic data.
# Run with: Rscript dev/issue-checks/issue-54.R
cat("Issue #54: ellipsoid_selection() fails with nvarstest = 1\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

set.seed(54)
env_train <- data.frame(bio_01 = rnorm(60, 20, 2), bio_12 = rnorm(60, 1000, 100))
env_test <- data.frame(bio_01 = rnorm(20, 21, 3), bio_12 = rnorm(20, 1050, 150))
env_bg <- data.frame(bio_01 = runif(1000, 5, 45), bio_12 = runif(1000, 200, 2500))
run <- function(nv) {
  tryCatch({
    out <- NULL
    utils::capture.output(
      out <- ellipsoid_selection(env_train = env_train, env_test = env_test,
                                 env_vars = c("bio_01", "bio_12"),
                                 nvarstest = nv, level = 0.975, mve = FALSE,
                                 env_bg = env_bg, omr_criteria = 0.5,
                                 parallel = FALSE, proc = FALSE))
    out
  }, error = function(e) paste("error:", conditionMessage(e)))
}
cat("EXPECTED: nvarstest = 1 returns a table with 2 one-variable models\n")
r1 <- run(1)
if (is.character(r1)) {
  cat("GOT:", r1, "\n")
  ok <- FALSE
} else {
  cat("GOT: table with", nrow(r1), "models; fitted_vars =",
      paste(r1$fitted_vars, collapse = " | "), "\n")
  ok <- nrow(r1) == 2 && all(r1$nvars == 1)
}
# Keep existing calls working: two variables still give one model
r2 <- run(2)
ok2 <- !is.character(r2) && nrow(r2) == 1
cat("Two-variable call still works:", ok2, "\n")
cat("\nRESULT:", if (ok && ok2) "PASS" else "FAIL", "\n")
