#' Function to thin longitude and latitude data
#' @description
#' Cleans up duplicated or redundant occurrence records that present overlapping
#' longitude and latitude coordinates. Thinning can be performed using either a
#' distance threshold (in coordinate units) or a pixel neighborhood approach.
#' @param data A data.frame with longitude and latitude of occurrence records.
#' @param longitude A character vector indicating the column name of the
#' "longitude" variable.
#' @param latitude A character vector indicating the column name of the
#' "latitude" variable.
#' @param threshold A numeric value representing the distance threshold between
#' coordinates to be considered duplicates, in the units of the coordinates,
#' which are the units of the layers' CRS (decimal degrees for lon/lat layers,
#' usually metres for projected layers). It is not converted to metres: for
#' lon/lat data, a threshold in degrees covers a shorter east-west distance at
#' higher latitudes. A common choice is the layer resolution, e.g.
#' \code{terra::res(raster_mask)[1]}. Used only when `by_mask` is
#' \code{FALSE}; with \code{by_mask = TRUE} the neighbourhood is set in
#' pixels by `n_ngbs`. Coordinates must be in the CRS of `raster_mask` (or of
#' the environmental layers); they are not reprojected.
#' @param by_mask Logical. If \code{T}, the thinning process will use a raster layer
#' as a mask for defining distance in pixel units.
#' @param raster_mask An object of class SpatRaster that serves as a reference
#' to thin the occurrence data. Required if `by_mask` is \code{T}.
#' @param n_ngbs Number of pixels used to define the neighborhood matrix that
#' helps determine which occurrences are duplicates:
#'   - 0 removes occurrences within the same pixel, keeping one.
#'   - 1 considers duplicates all occurrences within a distance of one pixel.
#'   - n considers duplicates all occurrences within a distance of n pixels.
#' @param process_ngbs_by Numeric parameter to improve memory management.
#' It process neighbor cells by a quantity specified by the user.
#' @param progress Logical. If \code{TRUE}, show computation progress.
#' @return Returns a data.frame with cleaned occurrence records, excluding
#' duplicates based on the specified criteria.
#' @details
#' This function cleans up duplicated occurrences based on the specified
#' distance threshold. If `by_mask` is \code{T}, the distance is interpreted as
#' pixel distance using the provided raster_mask; otherwise, it is interpreted
#' as a planar distance in coordinate units (degrees for lon/lat data). With
#' `by_mask = TRUE`, records outside the extent
#' of `raster_mask` are removed with a warning. Records with a missing
#' longitude or latitude are removed with a message.
#'
#' Distance thinning is greedy: records are visited from the last row to the
#' first, and a record is kept unless it lies within `threshold` of a record
#' that has already been kept. No two kept records are closer than
#' `threshold`, and every removed record is within `threshold` of a kept one.
#' @examples
#' data(abronia)
#' tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
#' tenm_mask <- terra::rast(file.path(tempora_layers_dir,"1939/bio_01.tif"))
#' # Clean duplicates without raster mask (just by distance threshold)
#' # First check the number of occurrence records
#' print(nrow(abronia))
#' # Clean duplicated records using a distance of ~ 18 km (0.1666667 grades)
#' ab_1 <- tenm::clean_dup(data =abronia,
#'                         longitude = "decimalLongitude",
#'                         latitude = "decimalLatitude",
#'                         threshold = terra::res(tenm_mask),
#'                         by_mask = FALSE,
#'                         raster_mask = NULL)
#' # Check number of records
#' print(nrow(ab_1))
#' # Clean duplicates using a raster mask
#' ab_2 <- tenm::clean_dup(data =abronia,
#'                         longitude = "decimalLongitude",
#'                         latitude = "decimalLatitude",
#'                         threshold = terra::res(tenm_mask)[1],
#'                         by_mask = TRUE,
#'                         raster_mask = tenm_mask,
#'                         n_ngbs = 1)
#' # Check number of records
#' print(nrow(ab_2))
#' @export
#'
clean_dup <- function(data,longitude,latitude,threshold=0.0, by_mask = FALSE,
                      raster_mask = NULL, n_ngbs = 0,
                      process_ngbs_by =100,progress =TRUE){
  if(isTRUE(by_mask) && !inherits(raster_mask, "SpatRaster"))
    stop("'raster_mask' must be a SpatRaster when by_mask = TRUE", call. = FALSE)
  na_xy <- is.na(data[[longitude]]) | is.na(data[[latitude]])
  if(any(na_xy))
    message(sum(na_xy), " record(s) with a missing longitude or latitude ",
            "were removed")
  data <- data[!na_xy, , drop = FALSE]

  #dat_sp <- sp::SpatialPointsDataFrame(data[,c(longitude ,latitude)],data)
  if(isTRUE(by_mask)){
    nbase <- 2*n_ngbs+1
    ngMat <- base::matrix(rep(1,nbase*nbase),
                          ncol =nbase,byrow = T )
    ngMat[n_ngbs+1,n_ngbs+1] <- 0
    cellids <- terra::cellFromXY(raster_mask,
                                 as.matrix(data[,c(longitude,latitude)]))
    # Records outside raster_mask have no cell; drop them in both branches
    n_out <- sum(is.na(cellids))
    if(n_out > 0L)
      warning(n_out, " record(s) fall outside the extent of 'raster_mask' and",
              " were removed. Check that the coordinates and 'raster_mask'",
              " use the same CRS.", call. = FALSE)
    ids_nodup <- which(!duplicated(cellids) & !is.na(cellids))
    cellids2 <- cellids[ids_nodup]
    # Subset the input so its columns and their order are kept
    dat2 <- data[ids_nodup, , drop = FALSE]
    if (n_ngbs == 0){
      return(dat2)
    } else {
      dat2_cells <- cellids2
      cellids2 <- base::sort(cellids2)
      n_cells <- length(cellids2)
      cut_offs <- ceiling(n_cells/process_ngbs_by)
      if(cut_offs > 1){
        pro_cells <- base::cut(x = seq_len(n_cells),cut_offs)
        cell_idsL <- base::split(cellids2,pro_cells)

      } else{
        cell_idsL <- list(cellids2)

      }


      #adj_cells <- terra::adjacent(x = raster_mask,cells=cellids2,
      #                             directions = ngMat,
      #                             pairs = TRUE)

      adj_cells <- seq_along(cell_idsL) |> purrr::map(function(x){
        adj_cells <- terra::adjacent(x = raster_mask,
                                     cells = cell_idsL[[x]],pairs=TRUE,
                                     directions = ngMat,include=FALSE)
        adj_cells <- data.frame(adj_cells)
        return(adj_cells)

      },.progress = progress)
      adj_cells <- purrr::compact(adj_cells)
      adj_cells <- purrr::list_rbind(adj_cells)

      #occ_adj_id <- which(adj_cells[,2] %in% cellids2)
      #if(length(occ_adj_id)>0L){
      #  adj_cells <- adj_cells[occ_adj_id,]
      #}

      adj_cellsL <- split(adj_cells[,2], adj_cells[,1])
      targets <- names(adj_cellsL)
      keep <- rep(NA,length(targets ))
      j <- 1
      for(i in seq_along(targets)){
        focal <- targets[i]
        if(focal %in% names(adj_cellsL)){
          vecinos <- as.character(adj_cellsL [[focal]])
          id_bas <- which(names(adj_cellsL) %in% vecinos)
          if(length(id_bas)>0L) adj_cellsL <- adj_cellsL[- id_bas]
          keep[j] <- as.numeric(focal)
          j <- j+1
        }

      }
      ids_nodup2 <- which(dat2_cells %in% keep)
      dat3 <- dat2[ids_nodup2, names(data), drop = FALSE]
      return(dat3)
    }

  }
  else{
    dat_sp <- sf::st_as_sf(data,coords=c(longitude,latitude))
    dat_bf <- sf::st_buffer(dat_sp, dist = threshold[1])
    # nbrs[[i]]: records whose threshold buffer contains record i
    nbrs <- sf::st_intersects(dat_sp,dat_bf)
    # Greedy thinning from the last record to the first: a record is kept
    # unless it is within the threshold of a record that is already kept.
    keep <- rep(TRUE, nrow(data))
    for(i in rev(seq_len(nrow(data)))){
      if(keep[i]){
        nb <- nbrs[[i]]
        keep[nb[nb < i]] <- FALSE
      }
    }
    no_dups <- which(keep)
    return(data[no_dups,])
  }
}
