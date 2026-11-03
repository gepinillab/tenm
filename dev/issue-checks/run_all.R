# Run every dev/issue-checks/issue-NN.R script and print a PASS/FAIL table.
#
# Usage (from the package root):
#   Rscript dev/issue-checks/run_all.R            # uses the installed tenm
#   Rscript dev/issue-checks/run_all.R <lib_dir>  # uses tenm installed in <lib_dir>
#
# Each script runs in its own R process. A script that stops with an error,
# or prints no "RESULT:" line, is reported as ERROR.
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else NULL

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else
  file.path(getwd(), "dev", "issue-checks")
scripts <- sort(list.files(here, pattern = "^issue-[0-9]+\\.R$", full.names = TRUE))
if (!length(scripts)) stop("No issue-NN.R scripts found in ", here)

rscript <- file.path(R.home("bin"), "Rscript")
env <- character()
if (!is.null(lib)) {
  old <- Sys.getenv("R_LIBS")
  env <- paste0("R_LIBS=", paste(c(lib, if (nzchar(old)) old), collapse = .Platform$path.sep))
}

rows <- lapply(scripts, function(s) {
  t0 <- Sys.time()
  out <- suppressWarnings(system2(rscript, shQuote(s), stdout = TRUE, stderr = TRUE, env = env))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  res <- grep("^RESULT:", out, value = TRUE)
  result <- if (length(res)) trimws(sub("^RESULT:", "", res[length(res)])) else "ERROR"
  ver <- grep("^tenm [0-9.]+ loaded from", out, value = TRUE)
  ver <- if (length(ver)) sub("^tenm ([0-9.]+) loaded from.*$", "\\1", ver[1]) else NA_character_
  title <- out[grepl("^Issue #[0-9]+:", out)]
  title <- if (length(title)) sub("^Issue #[0-9]+: *", "", title[1]) else ""
  if (nchar(title) > 70) title <- paste0(substr(title, 1, 67), "...")
  if (result == "ERROR") {
    cat("---- ", basename(s), " produced no RESULT line; last output:\n", sep = "")
    cat(utils::tail(out, 15), sep = "\n")
  }
  data.frame(issue = paste0("#", as.integer(sub("^issue-([0-9]+)\\.R$", "\\1", basename(s)))),
             result = result, tenm = ver, seconds = round(secs, 1), title = title)
})
tab <- do.call(rbind, rows)
options(width = 200)
cat("\n")
print(tab, row.names = FALSE, right = FALSE)
cat("\n", sum(tab$result == "PASS"), " PASS, ", sum(tab$result == "FAIL"), " FAIL, ",
    sum(!tab$result %in% c("PASS", "FAIL")), " ERROR\n", sep = "")
invisible(tab)
