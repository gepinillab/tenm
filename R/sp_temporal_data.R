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
#' @param occs A data.frame with information about occurrence records of the
#' species being modeled.
#' @param longitude Name of the variable in 'occs' containing longitude data.
#' @param latitude Name of the variable in 'occs' containing latitude data.
#' @param sp_date_var Name of the date variable.
#' @param occ_date_format Format of dates in occurrence records.
#' Options: "y", "ym", "ymd", "mdy", "my", "dmy".
#' @param layers_date_format Format of dates in raster layers. Options:
#' "y", "ym", "ymd", "mdy", "my", "dmy".
#' @param layers_by_date_dir Directory containing folders organized by date
#' with raster layers of environmental information.
#' @param layers_ext Extension or path of each raster layer archive
#' (e.g., ".tif").
#' @return Returns a  sp.temporal.modeling object (list) with the coordinates
#' of each occurrences points, the years of observation and the path to the
#' temporal layers.
#' @details
#'  The format of dates for each layer can be organized in a particular pattern,
#'  for example year/month/day ("ymd"), year/month ("ym"), just year ("y") or
#'  some other arrangement like month/year ("my"), month/year/day ("myd"),
#'  day/month/year ("dmy").
#'
#'  Each record is matched to the folder with the nearest date. A message
#'  reports how many records have no folder of the same date, and a warning
#'  reports records outside the range of folder dates. Records whose date is
#'  missing or cannot be read are removed with a warning.
#'
#'  The date folders used must all contain the same layer files. A warning is
#'  given when the layers have no CRS, or when records fall outside the extent
#'  of the layers (usually because the coordinates and the layers use
#'  different CRSs; occurrence coordinates are not reprojected).
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
                             layers_by_date_dir,layers_ext="*.tif$"){
  classes_sp_data <- c("data.frame")
  if(class(occs) %in% classes_sp_data){

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
      date_formats <- c("y","ym","ymd","my","myd","dmy")
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

      dates_ids <- nearest_date_id(occ_dates = date_sps,
                                   layers_dates = layers_dates)
      # Report records without a layer of the same date
      date_gap <- abs(as.numeric(date_sps - layers_dates[dates_ids]))
      n_occ <- length(date_gap)
      n_exact <- sum(date_gap == 0)
      if(n_exact < n_occ){
        message(n_exact, " of ", n_occ, " records match a layer date exactly; ",
                n_occ - n_exact, " were matched to the nearest layer date ",
                "(largest gap: ", max(date_gap), " days)")
      }
      rng <- range(layers_dates, na.rm = TRUE)
      outside <- date_sps < rng[1] | date_sps > rng[2]
      if(any(outside))
        warning(sum(outside), " record(s) are outside the range of layer dates (",
                rng[1], " to ", rng[2], ") and were matched to the first or ",
                "last layer date", call. = FALSE)

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
