# Run an issue check script against two versions of tenm.
#
# Usage (from inside the git repository):
#   Rscript dev/issue-checks/compare.R <issue> [base_ref] [pr_ref]
#   Rscript dev/issue-checks/compare.R 5 main batch/1-core-wrong-results
#   Rscript dev/issue-checks/compare.R all main batch/1-core-wrong-results
#
# <issue> is an issue number (e.g. 5) or "all". base_ref defaults to "main"
# and pr_ref to "HEAD". Each ref is exported with `git archive` and installed
# into its own temporary library; the check script is then run once per
# library, so the installed tenm in your normal library is never touched.
# The check scripts are read from pr_ref (falling back to the working tree),
# because they usually do not exist on the base ref.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript dev/issue-checks/compare.R <issue|all> [base_ref] [pr_ref]")
issue <- args[1]
base_ref <- if (length(args) >= 2) args[2] else "main"
pr_ref <- if (length(args) >= 3) args[3] else "HEAD"

git <- function(...) {
  out <- suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) stop("git ", paste(c(...), collapse = " "), " failed:\n",
                                            paste(out, collapse = "\n"))
  out
}
root <- git("rev-parse", "--show-toplevel")
resolve <- function(ref) {
  for (r in c(ref, paste0("origin/", ref))) {
    sha <- tryCatch(git("-C", root, "rev-parse", "--verify", "--quiet",
                        paste0(r, "^{commit}")), error = function(e) character())
    if (length(sha)) return(sha[1])
  }
  stop("Cannot resolve git ref '", ref, "'")
}
base_sha <- resolve(base_ref)
pr_sha <- resolve(pr_ref)

work <- tempfile("tenm-compare-")
dir.create(work)
on.exit(unlink(work, recursive = TRUE), add = TRUE)

# Collect the check scripts from pr_ref (or the working tree).
check_dir <- file.path(work, "checks")
dir.create(check_dir)
listed <- tryCatch(git("-C", root, "ls-tree", "--name-only", pr_sha, "dev/issue-checks/"),
                   error = function(e) character())
listed <- basename(grep("issue-[0-9]+\\.R$", listed, value = TRUE))
if (length(listed)) {
  for (f in listed) {
    writeLines(git("-C", root, "show", paste0(pr_sha, ":dev/issue-checks/", f)),
               file.path(check_dir, f))
  }
} else {
  local <- list.files(file.path(root, "dev", "issue-checks"), pattern = "^issue-[0-9]+\\.R$",
                      full.names = TRUE)
  file.copy(local, check_dir)
}
if (identical(tolower(issue), "all")) {
  scripts <- sort(list.files(check_dir, pattern = "^issue-[0-9]+\\.R$", full.names = TRUE))
} else {
  scripts <- file.path(check_dir, sprintf("issue-%02d.R", as.integer(issue)))
  if (!file.exists(scripts)) stop("No check script ", basename(scripts), " in ", pr_ref,
                                  " or the working tree")
}

install_ref <- function(sha, label) {
  src <- file.path(work, paste0("src-", label))
  lib <- file.path(work, paste0("lib-", label))
  dir.create(src)
  dir.create(lib)
  tarball <- file.path(work, paste0(label, ".tar"))
  git("-C", root, "archive", "--format=tar", "-o", tarball, sha)
  utils::untar(tarball, exdir = src)
  cat("Installing ", label, " (", substr(sha, 1, 7), ") into ", lib, "\n", sep = "")
  out <- suppressWarnings(system2(file.path(R.home("bin"), "R"),
                                  c("CMD", "INSTALL", "--no-test-load", "--no-docs",
                                    "--no-multiarch", "-l", shQuote(lib), shQuote(src)),
                                  stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) {
    cat(out, sep = "\n")
    stop("Installing ", label, " failed")
  }
  lib
}
libs <- c(base = install_ref(base_sha, "base"), pr = install_ref(pr_sha, "pr"))

run_check <- function(script, lib) {
  old <- Sys.getenv("R_LIBS")
  env <- paste0("R_LIBS=", paste(c(lib, if (nzchar(old)) old), collapse = .Platform$path.sep))
  # Run from a scratch directory so stray plot files do not land in the repo.
  owd <- setwd(work)
  on.exit(setwd(owd))
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
                                  stdout = TRUE, stderr = TRUE, env = env))
  res <- grep("^RESULT:", out, value = TRUE)
  list(out = out, result = if (length(res)) trimws(sub("^RESULT:", "", res[length(res)])) else "ERROR")
}

summary_rows <- list()
for (s in scripts) {
  results <- c()
  for (label in names(libs)) {
    ref <- if (label == "base") base_ref else pr_ref
    sha <- if (label == "base") base_sha else pr_sha
    cat("\n==== ", basename(s), " on ", label, ": ", ref, " (", substr(sha, 1, 7), ") ====\n", sep = "")
    r <- run_check(s, libs[[label]])
    keep <- grepl("^(Issue #|tenm |EXPECTED:|GOT:|NOTE:|RESULT:)", r$out)
    cat(if (r$result == "ERROR") utils::tail(r$out, 20) else r$out[keep], sep = "\n")
    results[label] <- r$result
  }
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    issue = paste0("#", as.integer(sub("^issue-([0-9]+)\\.R$", "\\1", basename(s)))),
    base = results[["base"]], pr = results[["pr"]],
    fixed = results[["base"]] != "PASS" && results[["pr"]] == "PASS")
}
options(width = 200)
cat("\n==== Summary: base = ", base_ref, " (", substr(base_sha, 1, 7), "), pr = ", pr_ref,
    " (", substr(pr_sha, 1, 7), ") ====\n", sep = "")
print(do.call(rbind, summary_rows), row.names = FALSE, right = FALSE)
