# Issue #27: ellipsoid_omr(): background prevalence includes the training
# presences.
# Synthetic data only. Run with: Rscript dev/issue-checks/issue-27.R
cat("Issue #27: ellipsoid_omr(): background prevalence includes the training presences\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

set.seed(27)
env_train <- data.frame(bio_01 = rnorm(60, 20, 3), bio_12 = rnorm(60, 1000, 150))
env_test <- data.frame(bio_01 = rnorm(20, 20, 3), bio_12 = rnorm(20, 1000, 150))
env_bg <- data.frame(bio_01 = runif(300, 5, 40), bio_12 = runif(300, 200, 2500))

r <- ellipsoid_omr(env_data = env_train, env_test = env_test, env_bg = env_bg,
                   cf_level = 0.975, mve = FALSE)

# Reference: share of background points inside the same ellipsoid.
m <- cov_center(env_train, mve = FALSE, level = 0.975, vars = 1:2)
in_bg <- inEllipsoid(m$centroid, m$covariance, env_bg, level = 0.975)$in_Ellipsoid
in_tr <- inEllipsoid(m$centroid, m$covariance, env_train, level = 0.975)$in_Ellipsoid
in_te <- inEllipsoid(m$centroid, m$covariance, env_test, level = 0.975)$in_Ellipsoid
prev_bg <- mean(in_bg)
prev_bg_train <- mean(c(in_bg, in_tr))
p_ref <- stats::binom.test(sum(in_te), length(in_te), p = prev_bg, alternative = "greater")$p.value

ok <- isTRUE(all.equal(r$bg_prevalence, prev_bg)) && isTRUE(all.equal(r$pval_bin, p_ref))

cat("EXPECTED: bg_prevalence = share of the", nrow(env_bg), "background points inside the ellipsoid =",
    format(prev_bg, digits = 6), "and pval_bin uses it =", format(p_ref, digits = 6), "\n")
cat("GOT:      bg_prevalence =", format(r$bg_prevalence, digits = 6),
    "(background + training presences would give", format(prev_bg_train, digits = 6), ")",
    "| pval_bin =", format(r$pval_bin, digits = 6), "\n")
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
