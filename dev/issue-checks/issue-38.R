# Issue #38: Move rgl from Imports to Suggests.
# Checks the installed DESCRIPTION, then runs a child R process in which rgl
# cannot be loaded (a broken dummy 'rgl' is put first on the library path):
# tenm must load, 2-D plots must work, and 3-D plots must stop with a clear
# error naming rgl.
# Run with: Rscript dev/issue-checks/issue-38.R
cat("Issue #38: Move rgl from Imports to Suggests\n")
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

desc <- utils::packageDescription("tenm")
in_imports <- grepl("\\brgl\\b", if (is.null(desc$Imports)) "" else desc$Imports)
in_suggests <- grepl("\\brgl\\b", if (is.null(desc$Suggests)) "" else desc$Suggests)

fake <- file.path(tempdir(), "issue38_fake_lib")
unlink(fake, recursive = TRUE)
dir.create(file.path(fake, "rgl"), recursive = TRUE)
writeLines(c("Package: rgl", "Version: 99.0", "Title: broken stand-in"),
           file.path(fake, "rgl", "DESCRIPTION"))
child <- tempfile(fileext = ".R")
writeLines(c(
  sprintf(".libPaths(c(%s, .libPaths()))", deparse(fake)),
  "cat('rgl loadable:', requireNamespace('rgl', quietly = TRUE), '\\n')",
  "ok <- tryCatch({ suppressPackageStartupMessages(library(tenm)); TRUE }, error = function(e) { cat('LOAD ERROR:', conditionMessage(e), '\\n'); FALSE })",
  "if (ok) {",
  "  grDevices::pdf(NULL)",
  "  set.seed(1); x <- rnorm(50); y <- rnorm(50); z <- rnorm(50)",
  "  r2 <- tryCatch({ plot_ellipsoid(x, y, mve = FALSE); 'ok' }, error = function(e) conditionMessage(e))",
  "  r3 <- tryCatch({ plot_ellipsoid(x, y, z, mve = FALSE); 'ok' }, error = function(e) conditionMessage(e))",
  "  cat('2D:', r2, '\\n'); cat('3D:', r3, '\\n')",
  "}"), child)
out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), shQuote(child),
                                stdout = TRUE, stderr = TRUE))
loaded <- !any(grepl("^LOAD ERROR", out))
two_d <- any(grepl("^2D: ok", out))
three_d <- grep("^3D:", out, value = TRUE)
three_d_ok <- length(three_d) == 1 && grepl("rgl", three_d) && grepl("install", three_d)

cat("EXPECTED: rgl in Suggests, not Imports; without a loadable rgl, tenm loads, a 2-D",
    "plot_ellipsoid() works and a 3-D one stops with an error naming rgl and how to install it\n")
cat("GOT:      rgl in Imports =", in_imports, "| in Suggests =", in_suggests,
    "| child process:", paste(trimws(out), collapse = " / "), "\n")
cat("RESULT:", if (!in_imports && in_suggests && loaded && two_d && three_d_ok) "PASS" else "FAIL", "\n")
