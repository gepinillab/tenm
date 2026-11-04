#' Function to create a Species Temporal Data object (STD object).
#' @description
#' Creates an object of class sp.temporal.modeling that contains
#' a list with four attributes:
#'   - temporal_df: A data frame with the following columns:
#'     - Longitude: Longitude coordinates of occurrence records.
#'     - Latitude: Latitude coordinates of occurrence records.
#'     - Date: Date variable indicating when the species were observed.
#'     - Layer Dates: Format of dates for each layer of environmental data.
#'     - Layers Path: Path to the bioclimatic layer corresponding to each year.
#'   - sp_date_var: Name of the date variable column in the occurrence records.
#'   - lon_lat_vars: Names of the longitude and latitude columns.
#'   - layers_ext: Final extension format of the environmental information
#'     (e.g., ".tif").
#' @param occs A data.frame (or tibble) with information about occurrence records of the
#' species being modeled.
#' @param longitude Name of the variable in 'occs' containing longitude data
#' (the x coordinate, in the CRS of the layers; see Details).
#' @param latitude Name of the variable in 'occs' containing latitude data
#' (the y coordinate, in the CRS of the layers; see Details).
#' @param sp_date_var Name of the date variable.
#' @param occ_date_format Format of dates in occurrence records.
#' Options: "y", "ym", "ymd", "mdy", "my", "myd", "dmy".
#' @param layers_date_format Format of dates in raster layers. Options:
#' "y", "ym", "ymd", "mdy", "my", "myd", "dmy".
#' @param layers_by_date_dir Directory containing folders organized by date
#' with raster layers of environmental information.
#' @param layers_ext Extension or path of each raster layer archive
#' (e.g., ".tif").
#' @param on_unmatched What to do with records whose date, floored to the
#' unit of the layer folders, has no layer folder: \code{"drop"} (default)
#' removes them, \code{"nearest"} matches them to the nearest layer date, and
#' \code{"error"} stops. See Details.
#' @return Returns a  sp.temporal.modeling object (list) with the coordinates
#' of each occurrences points, the years of observation and the path to the
#' temporal layers.
#' @details
#'  The format of dates for each layer can be organized in a particular pattern,
#'  for example year/month/day ("ymd"), year/month ("ym"), just year ("y") or
#'  some other arrangement like month/year ("my"), month/year/day ("myd"),
#'  day/month/year ("dmy").
#'
#'  Each occurrence date is floored to the unit of the layer folders (the
#'  year for \code{"y"} folders, the month for \code{"ym"} and \code{"my"},
#'  the day otherwise; e.g. \code{"1975-06-15"} becomes 1975 for \code{"y"}
#'  folders) and matched exactly to a folder date. \code{on_unmatched} sets
#'  what happens to records with no folder for their date: \code{"drop"}
#'  (default) removes them, \code{"nearest"} matches them to the folder with
#'  the nearest date (the behaviour of earlier versions) and \code{"error"}
#'  stops. A message reports how many records were dropped or moved, and
#'  their dates (and, for \code{"nearest"}, the folder date each was moved
#'  to); with \code{"nearest"}, a warning also counts records outside the
#'  range of folder dates. Records whose date is missing or cannot be read are
#'  removed with a warning. Because each record keeps only the layer folder of
#'  its own date, \code{\link[tenm]{clean_dup_by_date}} never compares records
#'  from different dates (unless \code{on_unmatched = "nearest"} moved them to
#'  the same folder).
#'
#'  Occurrence coordinates must be in the coordinate reference system (CRS)
#'  of the layers; they are never reprojected. For lon/lat layers (e.g.
#'  EPSG:4326) give longitude and latitude in decimal degrees; for projected
#'  layers, give x and y in the layers' units. The date folders used must all
#'  contain the same layer files, with the same extent, resolution and CRS;
#'  otherwise \code{sp_temporal_data()} stops. A warning is given when the
#'  layers have no CRS, or when records fall outside the extent of the layers
#'  (usually because the coordinates are not in the layers' CRS).
#'
#' @examples
#' library(tenm)
#' #A data.frame with occurrences points information of Abronia graminea.
#' # See help(abronia)
#' data("abronia")
#' tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
#' abt <- tenm::sp_temporal_data(occs = abronia,
#'                               longitude = "decimalLongitude",
#'                               latitude = "decimalLatitude",
#'                               sp_date_var = "year",
#'                               occ_date_format="y",
#'                               layers_date_format= "y",
#'                               layers_by_date_dir = tempora_layers_dir,
#'                               layers_ext="*.tif$")

#' @export

sp_temporal_data <- function(occs,longitude,
                             latitude,sp_date_var,
                             occ_date_format="y",
                             layers_date_format= "y",
                             layers_by_date_dir,layers_ext="*.tif$",
                             on_unmatched = c("drop", "nearest", "error")){
  on_unmatched <- match.arg(on_unmatched)
  if(inherits(occs, "data.frame")){
    # tibbles and data.tables are used as plain data.frames
    occs <- as.data.frame(occs)

    if(methods::is(occs, "data.frame")) {
      lon_lat_vars <- c(longitude,latitude)
      if(!all(lon_lat_vars %in% names(occs)))
        stop("\n Please provide species a valid 'longitude' and 'latitude'")
    }
    if(sp_date_var %in% names(occs)){

      nearest_date_id <- function(occ_dates,layers_dates){
        r1 <-  seq_along(occ_dates) |> purrr::map_int(function(x){
          id_layer_dir <- which.min(abs(occ_dates[x] - layers_dates))
          return(id_layer_dir)
        }) #|> do.call('c', .)
        return(r1)
      }
      date_formats <- c("y","ym","ymd","my","myd","dmy","mdy")
      if(!occ_date_format %in% date_formats)
        stop("occ_date_format should have one of these values ",
             paste(date_formats,collapse = " "))
      if(occ_date_format == date_formats[1]){
        date_in_occs <- lubridate::ymd(occs[,sp_date_var], truncated = 2L)
      } else if(occ_date_format == date_formats[2]){
        date_in_occs <- lubridate::ym(occs[,sp_date_var])
      } else if(occ_date_format == date_formats[3]){
        date_in_occs <- lubridate::ymd(occs[,sp_date_var])
      } else if(occ_date_format == date_formats[4]){
        date_in_occs <- lubridate::my(occs[,sp_date_var])
      } else if(occ_date_format == date_formats[5]){
        date_in_occs <- lubridate::myd(occs[,sp_date_var])
      } else if(occ_date_format == date_formats[6]){
        date_in_occs <- lubridate::dmy(occs[,sp_date_var])
      } else if(occ_date_format == date_formats[7]){
        date_in_occs <- lubridate::mdy(occs[,sp_date_var])
      }

      if(!dir.exists(layers_by_date_dir))
        stop("'layers_by_date_dir' does not exist: ", layers_by_date_dir,
             call. = FALSE)
      layers_names <- list.dirs(path = layers_by_date_dir,recursive = F,
                                full.names = F)
      if(!length(layers_names))
        stop("'layers_by_date_dir' has no date folders: ", layers_by_date_dir,
             call. = FALSE)

      if(!layers_date_format %in% date_formats)
        stop("layers_date_format should have one of these values ",
             paste(date_formats,collapse = " "))
      if(layers_date_format == date_formats[1]){
        layers_dates <- lubridate::ymd(layers_names, truncated = 2L)
      } else if(layers_date_format == date_formats[2]){
        layers_dates <- lubridate::ym(layers_names)
      } else if(layers_date_format == date_formats[3]){
        layers_dates <- lubridate::ymd(layers_names)
      } else if(layers_date_format == date_formats[4]){
        layers_dates <- lubridate::my(layers_names)
      } else if(layers_date_format == date_formats[5]){
        layers_dates <- lubridate::myd(layers_names)
      } else if(layers_date_format == date_formats[6]){
        layers_dates <- lubridate::dmy(layers_names)
      } else if(layers_date_format == date_formats[7]){
        layers_dates <- lubridate::mdy(layers_names)
      }

      if(all(is.na(layers_dates)))
        stop("No folder name in 'layers_by_date_dir' could be read with ",
             "layers_date_format = '", layers_date_format, "'. Folder names: ",
             paste(utils::head(layers_names, 5), collapse = ", "),
             call. = FALSE)
      if(anyNA(layers_dates))
        warning(sum(is.na(layers_dates)), " folder(s) in 'layers_by_date_dir' ",
                "could not be read with layers_date_format = '",
                layers_date_format, "' and were ignored: ",
                paste(utils::head(layers_names[is.na(layers_dates)], 5),
                      collapse = ", "), call. = FALSE)

      date_in_occs_index <- which(!is.na(date_in_occs))
      if(!length(date_in_occs_index))
        stop("No occurrence date in '", sp_date_var, "' could be read with ",
             "occ_date_format = '", occ_date_format, "'", call. = FALSE)
      if(length(date_in_occs_index) < length(date_in_occs))
        warning(length(date_in_occs) - length(date_in_occs_index),
                " record(s) with a missing date or a date that could not be ",
                "read with occ_date_format = '", occ_date_format,
                "' were removed", call. = FALSE)
      date_sps <- date_in_occs[date_in_occs_index]
      # environmental layers by date

      # Floor each occurrence date to the unit of the layer folders (e.g. a
      # "ymd" date to its year for "y" folders) and match it exactly.
      unit <- switch(layers_date_format, y = "year", ym = , my = "month", "day")
      unit_fmt <- switch(unit, year = "%Y", month = "%Y-%m", "%Y-%m-%d")
      date_floor <- lubridate::floor_date(date_sps, unit)
      dates_ids <- match(date_floor, layers_dates)
      unmatched <- is.na(dates_ids)
      n_occ <- length(dates_ids)
      if(any(unmatched)){
        tab_un <- table(format(date_floor[unmatched], unit_fmt))
        un_txt <- paste0(names(tab_un), " (", as.integer(tab_un), ")",
                         collapse = ", ")
        head_txt <- paste0(sum(unmatched), " of ", n_occ, " record(s) have no ",
                           "layer folder for their ", unit, ": ", un_txt)
        if(on_unmatched == "error"){
          stop(head_txt, ". Use on_unmatched = \"drop\" or \"nearest\"",
               call. = FALSE)
        } else if(on_unmatched == "drop"){
          if(all(unmatched))
            stop(head_txt, ". No record is left; use on_unmatched = ",
                 "\"nearest\" to match records to the nearest layer date",
                 call. = FALSE)
          message(head_txt, ". They were dropped (on_unmatched = \"drop\")")
        } else {
          # on_unmatched = "nearest": the nearest folder by date
          near <- nearest_date_id(occ_dates = date_sps[unmatched],
                                  layers_dates = layers_dates)
          dates_ids[unmatched] <- near
          moved <- paste(format(date_floor[unmatched], unit_fmt), "->",
                         format(layers_dates[near], unit_fmt))
          tab_mv <- table(moved)
          message(n_occ - sum(unmatched), " of ", n_occ, " records match a ",
                  "layer date exactly; ", sum(unmatched), " were matched to ",
                  "the nearest layer date (on_unmatched = \"nearest\"): ",
                  paste0(names(tab_mv), " (", as.integer(tab_mv), ")",
                         collapse = ", "))
          # A folder named by year (or month) covers that whole year (or month)
          rng <- range(layers_dates, na.rm = TRUE)
          last_end <- switch(layers_date_format,
                             y = rng[2] + lubridate::years(1),
                             ym = , my = rng[2] + months(1),
                             rng[2] + 1) - 1
          outside <- date_sps < rng[1] | date_sps > last_end
          if(any(outside))
            warning(sum(outside), " record(s) are outside the range of layer ",
                    "dates (", rng[1], " to ", last_end, ") and were matched ",
                    "to the first or last layer date", call. = FALSE)
        }
      }
      keep <- !is.na(dates_ids)
      date_in_occs_index <- date_in_occs_index[keep]
      dates_ids <- dates_ids[keep]

      layers_all <- list.dirs(path = layers_by_date_dir,recursive = FALSE,
                              full.names = TRUE)
      layers_all <- base::normalizePath(layers_all,winslash = "/")
      temporal_df <- data.frame(occs[date_in_occs_index,c(lon_lat_vars,sp_date_var)],
                                layer_dates = layers_dates[dates_ids],
                                layers_path=layers_all[dates_ids])
      validate_layers(temporal_df$layers_path, layers_ext,
                      temporal_df[, lon_lat_vars])

      sp_temp_data <- list(temporal_df = temporal_df,
                           sp_date_var = sp_date_var,
                           lon_lat_vars =lon_lat_vars ,
                           #layers_path = layers_all[dates_ids],
                           layers_ext = layers_ext
                           #multiband = multiband
                           )

      class(sp_temp_data) <- c("sp.temporal.modeling")

      return(sp_temp_data)
    }
    else
      stop("\n Please provide species a valid 'date_var'")
  }
  else
    stop("\n Please provide species occurrence data as a 'data.frame'")

}
