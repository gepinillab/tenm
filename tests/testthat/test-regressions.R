# Regression tests for confirmed bugs (see gepinillab/tenm issues).

testthat::test_that("#2: ellipsoid_selection() honours mve and proc_iter", {
  set.seed(42)
  env_train <- data.frame(bio_01 = c(rnorm(50, 20, 2), rnorm(10, 30, 6)),
                          bio_12 = c(rnorm(50, 1000, 100), rnorm(10, 1500, 300)))
  env_test <- data.frame(bio_01 = rnorm(20, 21, 3),
                         bio_12 = rnorm(20, 1050, 150))
  env_bg <- data.frame(bio_01 = runif(3000, 5, 45),
                       bio_12 = runif(3000, 200, 2500))
  set.seed(1)
  direct <- tenm::ellipsoid_omr(env_data = env_train, env_test = env_test,
                                env_bg = env_bg, cf_level = 0.975,
                                mve = FALSE, proc = TRUE, proc_iter = 7,
                                rseed = TRUE)
  set.seed(1)
  utils::capture.output(
    sel <- tenm::ellipsoid_selection(env_train = env_train,
                                     env_test = env_test,
                                     env_vars = c("bio_01", "bio_12"),
                                     nvarstest = 2, level = 0.975,
                                     mve = FALSE, env_bg = env_bg,
                                     omr_criteria = 0.1, parallel = FALSE,
                                     proc = TRUE, proc_iter = 7,
                                     rseed = TRUE))
  cols <- c("om_rate_train", "om_rate_test", "bg_prevalence", "pval_bin",
            "pval_proc", "env_bg_paucratio", "env_bg_auc")
  testthat::expect_equal(unlist(sel[1, cols]), unlist(direct[1, cols]))
})
