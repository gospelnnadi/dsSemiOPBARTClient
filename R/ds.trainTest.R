`%||%` <- function(a, b) if (is.null(a)) b else a

# dsSemiOPBARTTrainTest.R
# ---------------------------------------------------------------------------
# TRAIN (ds.semiOPBARTTrain(), an alias defined in dsSemiOPBARTClient.R) and
# TEST/PREDICT (this file) are separate functions, server- and client-side.
# Training never touches the held-out test object; prediction never re-fits
# or grows trees -- it only loads the FINAL global tree state via set_trees()
# and calls do_predict(), exactly mirroring predict_new_patients_semiopbart()'s
# documented approximation (f(x): point estimate from final trees; h(w) and
# thresholds: full retained posterior draws).
# ---------------------------------------------------------------------------

#source("ds.serialize.R")

# =====================  TEST / PREDICT  =====================================

# ---- client side ------------------------------------------------------------

# #' @param fit             output of ds.semiOPBARTTrain() -- must contain
# #'   $theta, $us, $data.name, $site_names
# #' @param data.name_test  the TEST object, distinct from the object `fit`
# #'   was trained on -- passing the same name is refused
# #' @export
# ds.semiOPBARTPredict <- function(fit, data.name_test = "semiOPBART_test",
#                                   newobj_pred = "semiOPBART_pred",
#                                   nfilter = 5, datasources = NULL, seed = 35) {
#   set.seed(seed)
#   if (identical(data.name_test, fit$data.name))
#     stop("data.name_test must differ from the object ds.semiOPBARTTrain() ",
#          "used -- refusing to evaluate on the training data")
#   if (is.null(datasources)) datasources <- DSI::datashield.connections_find()

#   untrained_requested <- setdiff(names(datasources), fit$site_names)
#   if (length(untrained_requested))
#     stop("ds.semiOPBARTPredict(): ", paste(untrained_requested, collapse = ", "),
#          " were not part of the ds.semiOPBARTTrain() run this `fit` came ",
#          "from. Under the current swap-based design every site keeps its ",
#          "own complete forest -- there is no single pooled 'global' tree ",
#          "set left to reconstruct and ship to a site that never trained. ",
#          "Drop these from `datasources`, or have them train too.")

#   DSI::datashield.aggregate(datasources,
#       call("semiOPBARTLocalPredictDS", semiOPBART_toSerialize(fit$theta),
#            semiOPBART_toSerialize(fit$us), data.name_test, newobj_pred, nfilter,
#            fit$state_name %||% ".semiOPBART_state", seed = seed))
# }


#' Predict at one or more sites -- MUST be sites that were part of the
#' `datasources` ds.semiOPBARTTrainE() actually ran across; there is no
#' never-trained-site case (see ds.semiOPBARTTrainE()'s docstring).
#' @param fit  output of ds.semiOPBARTTrainE()
#' @param prediction_method  "point" (default): use fit$theta_mean/us_mean
#'   -- one number per coefficient/threshold, cheaper to send and to
#'   compute with. "draws": use fit$theta_draws/us_draws in full --
#'   class probabilities are averaged over every retained posterior draw,
#'   same as Architecture D's own semiOPBARTLocalPredictDS() already does,
#'   giving properly-calibrated (typically wider) predictive uncertainty
#'   at the cost of a bigger payload and an O(n_draws) computation per
#'   prediction call. semiOPBARTLocalPredictEDS() (trainTestDS.R) branches
#'   on whether it receives a vector or a matrix -- no separate server
#'   function needed for each choice.
#' @export
ds.semiOPBARTPredictE <- function(fit, data.name_test,
                                  newobj_pred = "semiOPBART_pred_E",
                                  prediction_method = c( "draws","point"), # c("point", "draws"),
                                  nfilter = 5, datasources = NULL, seed = 35) {
  prediction_method <- match.arg(prediction_method)
  set.seed(seed)
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()

  untrained_requested <- setdiff(names(datasources), fit$site_names)
  if (length(untrained_requested))
    stop("ds.semiOPBARTPredictE(): ", paste(untrained_requested, collapse = ", "),
         " were not part of the ds.semiOPBARTTrainE() run this `fit` came ",
         "from -- Architecture E has no never-trained-site prediction path. ",
         "Drop these from `datasources`, or re-run ds.semiOPBARTTrainE() ",
         "including them.")
  if (identical(data.name_test, fit$data.name))
    stop("data.name_test must differ from the object ds.semiOPBARTTrainE() ",
         "used -- refusing to evaluate on the training data")

  theta_arg <- if (prediction_method == "draws") fit$theta_draws else fit$theta_mean
  us_arg    <- if (prediction_method == "draws") fit$us_draws    else fit$us_mean

  result <-DSI::datashield.aggregate(datasources,
      call("semiOPBARTLocalPredictEDS", semiOPBART_toSerialize(theta_arg),
           semiOPBART_toSerialize(us_arg), data.name_test, newobj_pred, nfilter,
           fit$state_name %||% ".semiOPBART_state", seed = seed))
     result 
}

#' Predict at one or more sites using Architecture F's model (from
#' ds.semiOPBARTCombineF()) -- point-estimate only. F never produces a
#' pooled DRAWS matrix (only site-level theta_mean/theta_cov get combined
#' analytically), so there is no "draws" choice to offer here the way
#' ds.semiOPBARTPredictE() can -- see ds_localCombine.R's header.
#' @param combined  output of ds.semiOPBARTCombineF()
#' @export
ds.semiOPBARTPredictF <- function(combined, data.name_test, outcome_col,
                                  newobj_pred = "semiOPBART_pred_F",
                                  nfilter = 5, datasources = NULL, seed = 35) {
  set.seed(seed)
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()

  trained <- names(combined$per_site_fits)
  untrained_requested <- setdiff(names(datasources), trained)
  if (length(untrained_requested))
    stop("ds.semiOPBARTPredictF(): ", paste(untrained_requested, collapse = ", "),
         " did not run ds.semiOPBARTFitF() -- Architecture F has no ",
         "never-trained-site prediction path. Drop these from `datasources`, ",
         "or have them run ds.semiOPBARTFitF() first.")

  result <- DSI::datashield.aggregate(datasources,
      call("semiOPBARTLocalPredictFDS", semiOPBART_toSerialize(combined$theta),
           semiOPBART_toSerialize(combined$us), data.name_test, outcome_col,
           newobj_pred, nfilter, combined$state_name %||% ".semiOPBART_local_fit_F", seed = seed))
result
}


