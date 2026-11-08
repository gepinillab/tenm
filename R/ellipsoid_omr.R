#' Compute omission rate and statistical metrics for ellipsoid models.
#'
#' @description
#' Computes omission rate and statistical metrics for ellipsoid models using
#' environmental data.
#' @param env_data A data frame containing the environmental data used for
#' modeling.
#' @param env_test A data frame with environmental testing data.
#' Default is NULL. If provided, the selection process includes p-values
#'  from a binomial test.
#' @param env_bg Environmental data sampled from the calibration area to compute
#'   the approximated prevalence of the model.
#' @param cf_level Proportion of points to be included in the ellipsoids.
#'   Equivalent to the error (E) proposed by Peterson et al. (2008).
#'   \doi{10.1016/j.ecolmodel.2007.11.008}.
#' @param mve Logical. If \code{TRUE}, computes a minimum volume ellipsoid using
#'   \code{\link[MASS]{cov.rob}} from the MASS package. If \code{FALSE},
#'   uses the covariance matrix of the input data.
#' @param proc Logical. If \code{TRUE}, performs a partial ROC test.
#' @param proc_iter Numeric. Total number of iterations for the partial ROC
#'   bootstrap.
#' @param rseed Logical. If \code{TRUE}, sets a random seed for the partial
#' ROC bootstrap. Default is \code{TRUE}.
#' @param E_percent Numeric value from 0 to 100 used as the threshold (E)
#' for the partial ROC test. Default is 5. See \code{\link[tenm]{pROC}}.
#' @param boost_percent Numeric value from 0 to 100 representing the
#' percentage of testing data used in each partial ROC bootstrap iteration.
#' Default is 50. See \code{\link[tenm]{pROC}}.
#' @param sub_sample Logical. Whether to use a subsample of the background
#' suitability values for the partial ROC test. Default is \code{TRUE}.
#' See \code{\link[tenm]{pROC}}.
#' @param sub_sample_size Size of the subsample used when \code{sub_sample}
#' is \code{TRUE}. Default is 1000. See \code{\link[tenm]{pROC}}.
#' @return A data.frame with the following columns:
#'   - "fitted_vars": Names of variables that were fitted.
#'   - "nvars": Number of fitted variables
#'   - "om_rate_train": Omission rate of the training data.
#'   - "non_pred_train_ids": Row IDs of non-predicted training data.
#'   - "om_rate_test"': Omission rate of the testing data.
#'   - "non_pred_test_ids": Row IDs of non-predicted testing data.
#'   - "bg_prevalence": Approximated prevalence of the model: the proportion
#'     of \code{env_bg} points inside the ellipsoid. It is also the null
#'     probability of the binomial test.
#'   - "pval_bin": p-value of the binomial test; \code{NA} if there are no
#'     testing records.
#'   - "pval_proc": p-value of the partial ROC test.
#'   - "env_bg_paucratio": Environmental partial AUC ratio value.
#'   - "env_bg_auc": Environmental AUC value.
#'
#'   If the ellipsoid cannot be fitted (for example, the covariance matrix is
#'   singular), a one-row data.frame with "fitted_vars", "nvars",
#'   "status" (\code{"fit_failed"}) and "fit_error" (the error message).

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
#' #This code is for running in parallel
#' future::plan("multisession",workers=2)
#' abex <- tenm::ex_by_date(this_species = abtc,train_prop=0.7)
#' abbg <- tenm::bg_by_date(this_species = abex,
#'                          buffer_ngbs=10,n_bg=50000)
#' future::plan("sequential")
#' edata <- abex$temporal_df
#' etrain <- edata[edata$trian_test=="Train",c("bio_05","bio_06","bio_12")]
#' etest <- edata[edata$trian_test=="Test",c("bio_05","bio_06","bio_12")]
#' bg <- abbg$env_bg[,c("bio_05","bio_06","bio_12")]
#' eor <- ellipsoid_omr(env_data=etrain,env_test=etest,env_bg=bg,
#'                      cf_level=0.975,proc=TRUE)
#' eor
#' }
#' @export

ellipsoid_omr <- function(env_data,env_test=NULL,env_bg,cf_level,mve=TRUE,
                          proc=FALSE,proc_iter=100,rseed=TRUE,
                          E_percent=5,boost_percent=50,
                          sub_sample=TRUE,sub_sample_size=1000){
  emd <- try(tenm::cov_center(data = env_data,
                               mve = mve,
                               level = cf_level,
                               vars = 1:ncol(env_data)),
             silent = TRUE)

  message1 <- attr(emd,"class")== "try-error"
  if(length(message1)>0L){
    # Report the failed fit instead of returning NULL
    return(data.frame(fitted_vars = paste(colnames(env_data), collapse = ","),
                      nvars = ncol(env_data),
                      status = "fit_failed",
                      fit_error = trimws(conditionMessage(attr(emd, "condition")))))
  }

  in_e <-  inEllipsoid(centroid = emd$centroid,
                       eShape = emd$covariance,
                       env_data = env_data,
                       level = cf_level)

  fails_train_ids <- which(in_e$in_Ellipsoid== 0)

  if(length(fails_train_ids)>0){
    fails_train_ids <- paste0(fails_train_ids,collapse = ",")
  } else {
    fails_train_ids <- NA
  }

  occs_table <- table( in_e$in_Ellipsoid)

  succsID <- which(names(occs_table) %in% "1")
  failsID <- which(names(occs_table) %in% "0")

  occs_succs <-  if(length(succsID)>0L){
    occs_table[[succsID]]
  } else{
    0
  }
  occs_fail <-  if(length(failsID)>0L){
    occs_table[[failsID]]
  } else{
    0
  }

  a_train <-  occs_fail
  omrate_train <- a_train /nrow( in_e)

  d_results <- data.frame(fitted_vars =paste(names(emd$centroid),
                                             collapse =  ","),
                          nvars=length(emd$centroid),
                          om_rate_train= omrate_train,
                          non_pred_train_ids = fails_train_ids)
  has_test <- is.data.frame(env_test) || is.matrix(env_test)
  if(has_test){
    in_etest <-  tenm::inEllipsoid(centroid = emd$centroid,
                                    eShape = emd$covariance,
                                    env_data = env_test,
                                    level = cf_level)

    fails_test_ids <- which(in_etest$in_Ellipsoid== 0)

    if(length(fails_test_ids)>0){
      fails_test_ids <- paste0(fails_test_ids,collapse = ",")
    } else {
      fails_test_ids <- NA
    }

    suits_val <- exp(-0.5*( in_etest$mh_dist))

    occs_table_test <- table(in_etest$in_Ellipsoid)

    succsID <- which(names(occs_table_test) %in% "1")
    failsID <- which(names(occs_table_test) %in% "0")

    occs_succs_test <-  if(length(succsID)>0L){
      occs_table_test[[succsID]]
    } else{
      0
    }
    occs_fail_test <-  if(length(failsID)>0L){
      occs_table_test[[failsID]]
    } else{
      0
    }
    a_test <-  occs_fail_test
    omrate_test <- a_test /nrow( in_etest)
    d_results <- data.frame(d_results,
                            om_rate_test=omrate_test,
                            non_pred_test_ids=fails_test_ids)
  }

  if(!is.null(env_bg)){

    env_bg <- data.frame(env_bg)
    in_ebg <-  tenm::inEllipsoid(centroid = emd$centroid,
                                  eShape = emd$covariance,
                                  env_data = env_bg,
                                  level = cf_level)
    suits_bg <- exp(-0.5*in_ebg$mh_dist)

    # Prevalence over the background only (not the training presences)
    bg_table <- table(in_ebg$in_Ellipsoid)
    succs_bg_ID <- which(names(bg_table) %in% "1")
    fails_bg_ID <- which(names(bg_table) %in% "0")

    bg_succs <-  if(length(succs_bg_ID)>0L){
      bg_table[[succs_bg_ID]]
    } else{
      0
    }

    bg_fails <-  if(length(fails_bg_ID)>0L){
      bg_table[[fails_bg_ID]]
    } else{
      0
    }
    prevBG <- bg_succs/(bg_fails+bg_succs)
    d_results <-data.frame( d_results,
                            bg_prevalence= prevBG)

    if(has_test){
      #bin_table <- table(c(in_ebg$in_Ellipsoid,
      #                     in_etest$in_Ellipsoid))
      #binBG <- bin_table[[2]]/(bin_table[[1]]+bin_table[[2]])
      test_fail <-  occs_fail_test
      test_succs <- occs_succs_test
      # One-sided binomial test: P(X >= test_succs)
      # With no test records there are no trials and no test: NA, not 1
      n_trials <- test_succs + test_fail
      p_bin <- if (n_trials == 0) {
        NA_real_
      } else {
        stats::pbinom(test_succs - 1,
                      size = n_trials,
                      prob = prevBG,
                      lower.tail = FALSE)
      }
      d_results <-data.frame( d_results,
                              pval_bin=p_bin)
      if(proc){
        proc <- try(tenm::pROC(suits_bg,test_data = suits_val,
                            n_iter = proc_iter,rseed = rseed,
                            E_percent = E_percent,
                            boost_percent = boost_percent,
                            sub_sample = sub_sample,
                            sub_sample_size = sub_sample_size))
        message1 <- attr(proc,"class")== "try-error"
        if(length(message1)>0L){
          pval_proc <- NA
          mean_aucratio <- NA
          mean_auc <- NA
        } else{
          pval_proc <- proc$pROC_summary[3]
          mean_aucratio <- proc$pROC_summary[2]
          mean_auc <- proc$pROC_summary[1]
        }

        d_results <-data.frame( d_results,
                                pval_proc,
                                env_bg_paucratio= mean_aucratio,
                                env_bg_auc = mean_auc)
      }


    }

  }
  # Keep the fitted ellipsoid so that callers can reuse the exact model
  # that was evaluated (see tenm_selection() and predict()).
  attr(d_results, "ellipsoid") <- list(centroid = emd$centroid,
                                       covariance = emd$covariance)
  return(d_results)
}
