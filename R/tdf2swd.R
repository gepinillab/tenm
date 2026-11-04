#' Temporal data.frame to Samples With Data format
#' @description
#' Converts a temporal data.frame to Samples With Data (SWD) table for use
#' with other modeling platforms such as MaxEnt.
#' @param this_species An object of class sp.temporal.env
#' (see \code{\link[tenm]{ex_by_date}} ) or sp.temporal.bg
#' (see \code{\link[tenm]{bg_by_date}}).
#' @param sp_name Character vector specifying the species name.
#' @return A data.frame formatted as Samples With Data (SWD) table. For an
#' sp.temporal.env object: \code{sp_name}, the longitude and latitude
#' columns, the date column and the environmental variables. For an
#' sp.temporal.bg object: \code{sp_name = "background"}, longitude,
#' latitude, \code{year} (the year of the layer date of each folder) and the
#' environmental variables. Columns are selected by name.
#' @export
#' @examples
#' \donttest{
#' library(tenm)
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
#' abtc <- tenm::clean_dup_by_date(abt,threshold = 10/60)
#' future::plan("multisession",workers=2)
#' abex <- tenm::ex_by_date(this_species = abtc,
#'                          train_prop=0.7)
#' abbg <- tenm::bg_by_date(abex,
#'                          buffer_ngbs=10,n_bg=50000)
#' future::plan("sequential")
#' # SWD table for occurrence records
#' occ_swd <- tdf2swd(this_species=abex,sp_name="abro_gram")
#' # SWD table for background data
#' bg_swd <- tdf2swd(this_species=abbg)
#' }

tdf2swd <- function(this_species,sp_name="sp"){
  cls <- class(this_species)
  if(length(cls)==1L && methods::is(this_species,"sp.temporal.env")){
    df_base <- this_species$temporal_df
    # Select columns by name: coordinates, date, environmental variables
    swd_df <- data.frame(sp_name,
                         df_base[, c(this_species$lon_lat_vars,
                                     this_species$sp_date_var,
                                     env_var_names(this_species)),
                                 drop = FALSE])
  } else if("sp.temporal.bg" %in% cls){
    df_base <- this_species$env_bg
    # Year from the layer date of each folder, not from the folder name
    tdf <- this_species$temporal_df
    years <- lubridate::year(tdf$layer_dates[match(df_base$ID_YEAR,
                                                   tdf$layers_path)])
    bg_vars <- setdiff(names(df_base),
                       c("ID_YEAR", this_species$lon_lat_vars))
    swd_df <- data.frame(sp_name="background",
                         df_base[, this_species$lon_lat_vars],
                         year=years,df_base[, bg_vars, drop = FALSE])
  } else{
    stop("this_species should be of class sp.temporal.env or sp.temporal.bg")
  }
  return(swd_df)

}

# Names of the environmental variables of an sp.temporal.env object. Objects
# created before "env_vars" was stored use the columns of env_data (without
# trian_test), or else the columns of temporal_df that are not metadata.
env_var_names <- function(this_species){
  if(!is.null(this_species$env_vars)) return(this_species$env_vars)
  if(!is.null(this_species$env_data))
    return(setdiff(names(this_species$env_data), "trian_test"))
  meta <- c(this_species$lon_lat_vars, this_species$sp_date_var,
            "layer_dates", "layers_path", "cell_ids_year", "trian_test")
  setdiff(names(this_species$temporal_df), meta)
}
