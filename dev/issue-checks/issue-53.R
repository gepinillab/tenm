# Issue #53: Use accessible, deterministic default palettes.
# plot_ellipsoid() and ellipsoid_projection(plot = TRUE) should not use the
# RNG, should default to viridis colours from grDevices::hcl.colors(), and
# should let users set the colours. Plots are drawn to temporary PNG files
# (2-D) or the rgl null device (3-D). Synthetic data.
# Run with: Rscript dev/issue-checks/issue-53.R
cat("Issue #53: Use accessible, deterministic default palettes\n")
options(rgl.useNULL = TRUE)
suppressPackageStartupMessages(library(tenm))
cat("tenm", format(utils::packageVersion("tenm")), "loaded from", find.package("tenm"), "\n\n")

set.seed(53)
x <- stats::rnorm(100)
y <- stats::rnorm(100)
z <- stats::rnorm(100)
png_of <- function(expr) {
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = 200, height = 200)
  r <- tryCatch({
    force(expr)
    TRUE
  }, error = function(e) conditionMessage(e))
  grDevices::dev.off()
  if (!isTRUE(r)) return(paste("error:", r))
  unname(tools::md5sum(f))
}
rng_kept <- function(expr) {
  set.seed(1)
  s0 <- .Random.seed
  r <- tryCatch({
    suppressWarnings(force(expr))
    TRUE
  }, error = function(e) FALSE)
  isTRUE(r) && identical(s0, .Random.seed)
}

# (a) plot_ellipsoid(): RNG untouched; same picture whatever the seed; the
# default is hcl.colors(1, "viridis"); col is honoured
a_rng <- rng_kept(png_of(plot_ellipsoid(x, y, mve = FALSE)))
set.seed(1)
p1 <- png_of(plot_ellipsoid(x, y, mve = FALSE))
set.seed(2)
p2 <- png_of(plot_ellipsoid(x, y, mve = FALSE))
p_vir <- png_of(plot_ellipsoid(x, y, mve = FALSE, col = grDevices::hcl.colors(1, "viridis")))
p_red <- png_of(plot_ellipsoid(x, y, mve = FALSE, col = "red"))
ok_a <- a_rng && identical(p1, p2) && identical(p1, p_vir) && !identical(p1, p_red)

# (b) ellipsoid_projection(plot = TRUE), 2 variables: default equals viridis
# and a user palette changes the plot
mk_layer <- function() terra::rast(nrows = 20, ncols = 20, vals = stats::runif(400, 1, 3))
r2 <- c(mk_layer(), mk_layer())
names(r2) <- c("a", "b")
cen2 <- c(a = 2, b = 2)
cov2 <- diag(c(0.2, 0.2), 2)
dimnames(cov2) <- list(names(cen2), names(cen2))
b_rng <- rng_kept(png_of(ellipsoid_projection(r2, cen2, cov2, plot = TRUE, size = 2)))
q_def <- png_of(ellipsoid_projection(r2, cen2, cov2, plot = TRUE, size = 2))
q_vir <- png_of(ellipsoid_projection(r2, cen2, cov2, plot = TRUE, size = 2,
                                     col = grDevices::hcl.colors(100, "viridis")))
q_grey <- png_of(ellipsoid_projection(r2, cen2, cov2, plot = TRUE, size = 2,
                                      col = grDevices::grey.colors(10)))
ok_b <- b_rng && identical(q_def, q_vir) && !identical(q_def, q_grey)

# (c) 3 variables (rgl null device): RNG untouched and col accepted
r3 <- c(r2, mk_layer())
names(r3) <- c("a", "b", "c")
cen3 <- c(a = 2, b = 2, c = 2)
cov3 <- diag(0.2, 3)
dimnames(cov3) <- list(names(cen3), names(cen3))
c_rng <- rng_kept(ellipsoid_projection(r3, cen3, cov3, plot = TRUE, size = 2))
c_col <- tryCatch({
  ellipsoid_projection(r3, cen3, cov3, plot = TRUE, size = 2, col = grDevices::grey.colors(10))
  "accepted"
}, error = function(e) paste("error:", conditionMessage(e)))
ok_c <- c_rng && identical(c_col, "accepted")

cat("EXPECTED: (a) plot_ellipsoid() leaves .Random.seed unchanged, draws the same picture",
    "after set.seed(1) and set.seed(2), its default equals col = hcl.colors(1, 'viridis'),",
    "and col = 'red' changes it; (b) ellipsoid_projection() 2-D plot leaves the RNG unchanged,",
    "its default equals col = hcl.colors(100, 'viridis') and a grey palette changes it;",
    "(c) 3-D plot leaves the RNG unchanged and accepts col\n")
cat("GOT:      (a) RNG unchanged =", a_rng, "| same across seeds =", identical(p1, p2),
    "| default == viridis =", identical(p1, p_vir), "| col honoured =", !identical(p1, p_red), "\n")
cat("          (b) RNG unchanged =", b_rng, "| default == viridis =", identical(q_def, q_vir),
    "| col honoured =", !identical(q_def, q_grey), "\n")
cat("          (c) RNG unchanged =", c_rng, "| col:", c_col, "\n")
cat("RESULT:", if (ok_a && ok_b && ok_c) "PASS" else "FAIL", "\n")
