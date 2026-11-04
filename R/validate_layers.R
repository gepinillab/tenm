# Internal: check the date folders used by an sp.temporal.modeling object.
# Stops when a folder has no layers or lacks a layer found in the first
# folder; warns when the layers have no CRS or when occurrences fall outside
# the layers' extent (usually a CRS mismatch). Stops when a layer's grid
# (extent, resolution, CRS) differs from the first layer. Occurrence
# coordinates are never reprojected: they must use the layers' CRS.
validate_layers <- function(layers_path, layers_ext, xy){
  folders <- unique(layers_path)
  files <- lapply(folders, function(f){
    list.files(f, pattern = layers_ext, full.names = TRUE)
  })
  no_files <- lengths(files) == 0L
  if(any(no_files))
    stop("No layers matching layers_ext = '", layers_ext, "' in: ",
         paste(folders[no_files], collapse = ", "), call. = FALSE)
  lnames <- lapply(files, function(x) sort(basename(x)))
  ref <- lnames[[1]]
  differ <- vapply(lnames, function(x) !identical(x, ref), logical(1))
  if(any(differ)){
    detail <- vapply(which(differ), function(i){
      missing <- setdiff(ref, lnames[[i]])
      extra <- setdiff(lnames[[i]], ref)
      paste0(basename(folders[i]), ":",
             if(length(missing)) paste0(" missing ", paste(missing, collapse = ", ")),
             if(length(extra)) paste0(" extra ", paste(extra, collapse = ", ")))
    }, character(1))
    stop("All date folders must contain the same layers as '",
         basename(folders[1]), "' (", paste(ref, collapse = ", "), "):\n  ",
         paste(utils::head(detail, 10), collapse = "\n  "), call. = FALSE)
  }
  r1 <- terra::rast(files[[1]][1])
  # Every layer must share the grid (extent, rows/columns, resolution and
  # CRS) of the first layer of the first folder.
  all_files <- unlist(files, use.names = FALSE)
  same_geom <- vapply(all_files[-1], function(x){
    isTRUE(tryCatch(terra::compareGeom(r1, terra::rast(x)),
                    error = function(e) FALSE))
  }, logical(1))
  if(!all(same_geom))
    stop("All layers in the date folders must share the extent, resolution ",
         "and CRS of '", all_files[1], "'. These differ:\n  ",
         paste(utils::head(all_files[-1][!same_geom], 10), collapse = "\n  "),
         call. = FALSE)
  no_crs <- vapply(files, function(x) !nzchar(terra::crs(terra::rast(x[1]))),
                   logical(1))
  if(any(no_crs))
    warning("The layers in ", sum(no_crs), " date folder(s) have no CRS; ",
            "occurrence coordinates are assumed to use the layers' coordinate ",
            "system", call. = FALSE)
  cells <- terra::cellFromXY(r1, as.matrix(xy))
  n_out <- sum(is.na(cells))
  if(n_out > 0L)
    warning(n_out, " of ", length(cells), " record(s) fall outside the extent ",
            "of the layers. Check that the occurrence coordinates use the ",
            "layers' CRS (", terra::crs(r1, describe = TRUE)$name, ")",
            call. = FALSE)
  invisible(TRUE)
}
