# Writes small synthetic year folders for regression tests. In folder i,
# bio_01 = i * 10000 + cell and bio_12 = i * 10000 + 5000 + cell, so any
# extracted value identifies the folder and the cell it was read from.
make_year_layers <- function(folders, nrows = 30, ncols = 30,
                             crop_folder = NULL) {
  root <- tempfile("tenm_layers_")
  dir.create(root)
  tmpl <- terra::rast(nrows = nrows, ncols = ncols, xmin = 0, xmax = ncols,
                      ymin = 0, ymax = nrows, crs = "EPSG:4326")
  for (i in seq_along(folders)) {
    dir.create(file.path(root, folders[i]))
    r1 <- tmpl
    terra::values(r1) <- i * 10000 + seq_len(terra::ncell(tmpl))
    r2 <- tmpl
    terra::values(r2) <- i * 10000 + 5000 + seq_len(terra::ncell(tmpl))
    if (identical(folders[i], crop_folder)) {
      e <- terra::ext(tmpl)
      r1 <- terra::crop(r1, terra::ext(e[1] + 3, e[2], e[3], e[4]))
      r2 <- terra::crop(r2, terra::ext(e[1] + 3, e[2], e[3], e[4]))
    }
    names(r1) <- "bio_01"
    names(r2) <- "bio_12"
    terra::writeRaster(r1, file.path(root, folders[i], "bio_01.tif"))
    terra::writeRaster(r2, file.path(root, folders[i], "bio_12.tif"))
  }
  list(root = root, template = tmpl)
}
