test_that("predict forwards size to ellipsoid_projection", {
  set.seed(1)
  temporal_df <- data.frame(
    x = stats::rnorm(30),
    y = stats::rnorm(30),
    trian_test = "Train"
  )
  object <- structure(
    list(
      temporal_df = temporal_df,
      env_bg = data.frame(x = 0, y = 0),
      mods_table = data.frame(fitted_vars = "x,y")
    ),
    class = "sp.temporal.selection"
  )
  layers <- terra::rast(
    nrows = 2, ncols = 2, nlyrs = 2,
    xmin = 0, xmax = 2, ymin = 0, ymax = 2
  )
  names(layers) <- c("x", "y")
  terra::values(layers) <- matrix(stats::rnorm(8), ncol = 2)

  raster <- predict(object, layers = layers, mve = FALSE, size = 3)
  expect_identical(class(raster)[1], "SpatRaster")

  raster_list <- predict(object, layers = list(layers), mve = FALSE, size = 3)
  expect_identical(class(raster_list), "list")
  expect_identical(class(raster_list[[1]][[1]])[1], "SpatRaster")

  layer_dir <- tempfile()
  dir.create(layer_dir)
  on.exit(unlink(layer_dir, recursive = TRUE), add = TRUE)
  terra::writeRaster(layers[[1]], file.path(layer_dir, "x.tif"))
  terra::writeRaster(layers[[2]], file.path(layer_dir, "y.tif"))
  path_raster <- predict(
    object, layers_path = layer_dir, layers_ext = "\\.tif$",
    mve = FALSE, size = 3
  )
  expect_identical(class(path_raster)[1], "SpatRaster")
})
