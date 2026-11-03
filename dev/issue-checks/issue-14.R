# Issue #14: ellipsoid_omr(): exists("in_etest") matches objects in the
# global environment.
# Synthetic data only. Run with: Rscript dev/issue-checks/issue-14.R
cat("Issue #14: ellipsoid_omr(): exists(\"in_etest\") matches objects in the global environment\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

set.seed(14)
env_train <- data.frame(bio_01 = rnorm(60, 20, 3), bio_12 = rnorm(60, 1000, 150))
env_bg <- data.frame(bio_01 = runif(2000, 5, 40), bio_12 = runif(2000, 200, 2500))

# Reference: no object called in_etest anywhere.
ref <- ellipsoid_omr(env_data = env_train, env_test = NULL, env_bg = env_bg,
                     cf_level = 0.975, mve = FALSE)

# Same call with an unrelated user object called in_etest in the global env.
in_etest <- "anything"
got <- tryCatch(ellipsoid_omr(env_data = env_train, env_test = NULL, env_bg = env_bg,
                              cf_level = 0.975, mve = FALSE),
                error = function(e) e)
rm(in_etest)

ok <- is.data.frame(got) && identical(names(got), names(ref)) &&
  isTRUE(all.equal(got, ref))

cat("EXPECTED: with env_test = NULL and a global object 'in_etest', ellipsoid_omr()",
    "returns the same table as without it (columns:", paste(names(ref), collapse = ", "), ")\n")
cat("GOT:     ", if (inherits(got, "error")) paste("error:", conditionMessage(got)) else
  paste("columns:", paste(names(got), collapse = ", ")), "\n")
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
