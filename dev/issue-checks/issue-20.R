# Issue #20: ellipsoid_omr(): non_pred_test_ids is "" instead of NA (tests
# the wrong variable).
# Synthetic data only. Run with: Rscript dev/issue-checks/issue-20.R
cat("Issue #20: ellipsoid_omr(): non_pred_test_ids is \"\" instead of NA (tests the wrong variable)\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

set.seed(20)
env_train <- data.frame(bio_01 = rnorm(60, 20, 3), bio_12 = rnorm(60, 1000, 150))
# Test points close to the centroid (all inside) and one far away (outside).
env_test_in <- data.frame(bio_01 = c(20, 20.5, 19.5, 21, 19), bio_12 = c(1000, 990, 1010, 1020, 980))
env_test_out <- rbind(env_test_in, data.frame(bio_01 = 60, bio_12 = 5000))
env_bg <- data.frame(bio_01 = runif(2000, 5, 40), bio_12 = runif(2000, 200, 2500))

r_in <- ellipsoid_omr(env_data = env_train, env_test = env_test_in, env_bg = env_bg,
                      cf_level = 0.975, mve = FALSE)
r_out <- ellipsoid_omr(env_data = env_train, env_test = env_test_out, env_bg = env_bg,
                       cf_level = 0.975, mve = FALSE)

ok <- r_in$om_rate_test == 0 && is.na(r_in$non_pred_test_ids) &&
  identical(as.character(r_out$non_pred_test_ids), "6")

cat("EXPECTED: 0 of 5 test points omitted -> non_pred_test_ids = NA;",
    "point 6 omitted -> non_pred_test_ids = \"6\"\n")
cat("GOT:      om_rate_test =", r_in$om_rate_test, "-> non_pred_test_ids =",
    if (is.na(r_in$non_pred_test_ids)) "NA" else dQuote(r_in$non_pred_test_ids, FALSE),
    "; om_rate_test =", format(r_out$om_rate_test, digits = 4), "-> non_pred_test_ids =",
    dQuote(as.character(r_out$non_pred_test_ids), FALSE), "\n")
cat("RESULT:", if (isTRUE(ok)) "PASS" else "FAIL", "\n")
