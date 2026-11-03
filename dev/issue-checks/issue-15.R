# Issue #15: tenm_selection(): RandomPercent, sub_sample and sub_sample_size
# are never used.
# Synthetic data only. Run with: Rscript dev/issue-checks/issue-15.R
cat("Issue #15: tenm_selection(): RandomPercent, sub_sample and sub_sample_size are never used\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

# A minimal sp.temporal.bg object: tenm_selection() only reads temporal_df
# (with the trian_test column) and env_bg.
set.seed(15)
tdf <- data.frame(bio_01 = rnorm(80, 20, 3), bio_12 = rnorm(80, 1000, 150),
                  trian_test = rep(c("Train", "Test"), c(60, 20)))
sp_bg <- list(temporal_df = tdf, sp_date_var = "year", lon_lat_vars = c("x", "y"),
              layers_ext = ".tif$",
              env_bg = data.frame(bio_01 = runif(3000, 5, 40),
                                  bio_12 = runif(3000, 200, 2500)))
class(sp_bg) <- "sp.temporal.bg"

# Record the arguments pROC() receives.
rec <- new.env()
rec$calls <- list()
invisible(suppressMessages(trace("pROC",
  tracer = quote(assign("calls", c(get("calls", envir = rec),
                                   list(list(boost_percent = boost_percent,
                                             sub_sample = sub_sample,
                                             sub_sample_size = sub_sample_size,
                                             E_percent = E_percent,
                                             n_iter = n_iter))),
                        envir = rec)),
  where = asNamespace("tenm"), print = FALSE)))

want <- list(boost_percent = 30, sub_sample = FALSE, sub_sample_size = 200,
             E_percent = 10, n_iter = 9)
args <- list(this_species = sp_bg, omr_criteria = 0.1, ellipsoid_level = 0.975,
             vars2fit = c("bio_01", "bio_12"), nvars_to_fit = 2, mve = FALSE,
             proc = TRUE, RandomPercent = 30, sub_sample = FALSE,
             sub_sample_size = 200, NoOfIteration = 9, parallel = FALSE)
has_E <- "E_percent" %in% names(formals(tenm_selection))
if (has_E) args$E_percent <- 10
invisible(utils::capture.output(sel <- do.call(tenm_selection, args)))
suppressMessages(untrace("pROC", where = asNamespace("tenm")))

got <- rec$calls[[1]]
fmt <- function(x) paste(names(x), vapply(x, format, ""), sep = " = ", collapse = ", ")
ok <- has_E && length(rec$calls) == 1 && isTRUE(all.equal(got[names(want)], want))

cat("EXPECTED: tenm_selection(RandomPercent = 30, sub_sample = FALSE, sub_sample_size = 200,",
    "E_percent = 10, NoOfIteration = 9) calls pROC() with", fmt(want), "\n")
cat("GOT:      tenm_selection()", if (has_E) "accepts" else "has no", "E_percent argument;",
    "pROC() received", fmt(got[names(want)]), "\n")
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
