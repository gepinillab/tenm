# Issue #13: pROC() ignores rseed and resets the global RNG state.
# Synthetic suitability values only. Run with: Rscript dev/issue-checks/issue-13.R
cat("Issue #13: pROC() ignores rseed and resets the global RNG state\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

set.seed(1)
bgv <- runif(5000)
tv <- runif(30, .3, 1)

# A: rseed = FALSE must not reseed: the caller's stream continues.
set.seed(111); invisible(pROC(bgv, tv, n_iter = 20, rseed = FALSE)); a <- runif(1)
set.seed(999); invisible(pROC(bgv, tv, n_iter = 20, rseed = FALSE)); b <- runif(1)
ok_a <- !identical(a, b)

# B: rseed = FALSE bootstraps come from the caller's stream, so two different
# caller seeds give different bootstrap results (sub_sample = FALSE, so the
# bootstrap is the only random step).
set.seed(111); r111 <- pROC(bgv, tv, n_iter = 20, rseed = FALSE, sub_sample = FALSE)
set.seed(999); r999 <- pROC(bgv, tv, n_iter = 20, rseed = FALSE, sub_sample = FALSE)
ok_b <- !identical(r111$pROC_results, r999$pROC_results)

# C: rseed = TRUE is reproducible and leaves the caller's RNG state unchanged.
set.seed(222); before <- .Random.seed
t1 <- pROC(bgv, tv, n_iter = 20, rseed = TRUE)
state_kept <- identical(.Random.seed, before)
set.seed(222); t2 <- pROC(bgv, tv, n_iter = 20, rseed = TRUE)
ok_c <- state_kept && identical(t1$pROC_results, t2$pROC_results)

cat("EXPECTED: A) runif(1) after pROC(rseed = FALSE) differs between set.seed(111) and",
    "set.seed(999); B) rseed = FALSE bootstrap results differ between those seeds;",
    "C) pROC(rseed = TRUE) is reproducible and leaves .Random.seed unchanged\n")
cat("GOT:      A) runif(1) =", format(a, digits = 7), "vs", format(b, digits = 7),
    "; B) bootstrap results", if (ok_b) "differ" else "identical",
    "; C) reproducible:", identical(t1$pROC_results, t2$pROC_results),
    "| .Random.seed unchanged:", state_kept, "\n")
cat("RESULT:", if (ok_a && ok_b && ok_c) "PASS" else "FAIL", "\n")
