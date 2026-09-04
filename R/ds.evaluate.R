# # dsSemiOPBARTEvaluate.R
# # ---------------------------------------------------------------------------
# # Federated evaluation, mirroring evaluate_semiopbart()'s full metric set
# # (accuracy, balanced accuracy, macro F1, kappa, precision, recall, plus the
# # binary 0-vs-non-0 versions and both confusion matrices).
# #
# # KEY DESIGN POINT: every one of those metrics is a function of the
# # confusion matrix alone -- and confusion-matrix cells are just counts,
# # which sum safely across sites (subject to the usual nfilter floor). So
# # unlike training, evaluation does NOT need an iterative round-trip
# # protocol at all: one local confusion-matrix computation per site, one
# # aggregation, done. This file is intentionally generic -- it works
# # identically whether `newobj_pred` was produced by
# # dsSemiOPBARTTrainTest.R (Architecture D) or by the predict addendum in
# # dsSemiOPBARTLocalCombine.R (Option C), which is what makes the two
# # directly comparable in the orchestrator.
# # ---------------------------------------------------------------------------

# # ---- client side ------------------------------------------------------------
# #' Ordinal mean absolute error computed DIRECTLY from an aggregated
# #' (truth x predicted) confusion matrix, i.e. from counts that are already
# #' safe to have crossed the DataSHIELD boundary (see ds.semiOPBARTEvaluate()
# #' / semiOPBARTLocalConfusionDS()) -- no new disclosure surface, and no
# #' reliance on raw per-patient predictions ever reaching the client (they
# #' never do, and never should).
# #'
# #' cm's rows and columns are both ordered by `levels` (see
# #' semiOPBARTLocalConfusionDS(): `factor(..., levels = levels)` for both
# #' truth and prediction), so `outer(lv, lv, "-")` lines up with cm exactly.
# .semiOPBART_ordinalMAE <- function(cm, levels) {
#   lv <- suppressWarnings(as.numeric(as.character(levels)))
#   if (anyNA(lv)) lv <- seq_along(levels) - 1  # fallback: treat as equally-spaced ranks
#   dist_mat <- abs(outer(lv, lv, "-"))
#   sum(cm * dist_mat) / sum(cm)
# }

# #' Aggregate confusion matrices across sites and compute the full
# #' evaluate_semiopbart()-equivalent metric set from the pooled matrix.
# #'
# #' @param pred_obj  name of the local prediction object at each site
# #' @param levels    full ordinal level set (same as used in prepare/predict)
# #' @param label     a label for this evaluation run (e.g. "Architecture D",
# #'                  "Option C"), carried through for the comparison step
# #' @param validation_col, validation_data.name  see
# #'   semiOPBARTLocalConfusionDS(). If validation_col is given, sites
# #'   without that column (or without enough non-missing values) are simply
# #'   excluded from the aggregate, with a message naming them, rather than
# #'   failing the whole evaluation.
# #' @param dichotomize_threshold  NULL (default): binary evaluation uses the
# #'   original baseline-vs-rest split. A number (e.g. 2): binary evaluation
# #'   dichotomizes at this threshold instead (truth/prediction >= threshold
# #'   -> "1"). See semiOPBARTLocalConfusionDS() for the full rationale.
# #' @export
# ds.semiOPBARTEvaluate <- function(pred_obj = "semiOPBART_pred", levels,
#                                    nfilter = 5, datasources = NULL,
#                                    label = "model",
#                                    validation_col = NULL,
#                                    validation_data.name = NULL,
#                                    dichotomize_threshold = NULL) {
#   if (is.null(datasources)) datasources <- DSI::datashield.connections_find()

#   site_cm_raw <- DSI::datashield.aggregate(datasources,
#       call("semiOPBARTLocalConfusionDS", pred_obj, semiOPBART_toSerialize(levels), nfilter,
#            if (is.null(validation_col)) "null" else validation_col,
#            if (is.null(validation_data.name)) "null" else validation_data.name,
#            if (is.null(dichotomize_threshold)) "null" else as.character(dichotomize_threshold)))

#   no_validation <- names(site_cm_raw)[vapply(site_cm_raw, is.null, logical(1))]
#   if (length(no_validation) && !is.null(validation_col))
#     message("ds.semiOPBARTEvaluate(): no '", validation_col, "' validation ",
#             "available at: ", paste(no_validation, collapse = ", "),
#             " -- excluded from this evaluation, not treated as an error.")

#   site_cm <- site_cm_raw[!vapply(site_cm_raw, is.null, logical(1))]
#   if (length(site_cm) == 0)
#     stop("no site could be validated -- ",
#          if (!is.null(validation_col))
#            paste0("'", validation_col, "' was unavailable everywhere")
#          else "no predictions with a training-label truth were found")

#   site_cm <- site_cm_raw
# ####################–––––––––––––––––––##########################
#   # Aggregation uses the UNMASKED $confusion (not $confusion_reported) so
#   # cross-site sums are exact even when an individual site's cell was too
#   # small to report on its own -- the pooled cell is what actually needs
#   # the nfilter check, not each site's contribution to it.
#   cm     <- Reduce(`+`, lapply(site_cm, `[[`, "confusion"))
#   cm_bin <- Reduce(`+`, lapply(site_cm, `[[`, "confusion_bin"))
#   if (any(cm > 0 & cm < nfilter))
#     warning("a pooled multiclass confusion cell is below nfilter=", nfilter,
#             " even after aggregation -- treat that cell's contribution to ",
#             "per-class metrics with caution")

#   n <- sum(cm)
#   acc <- sum(diag(cm)) / n
#   mae <- .semiOPBART_ordinalMAE(cm, levels)
#   recall_per_class <- diag(cm) / rowSums(cm)
#   precision_per_class <- diag(cm) / colSums(cm)
#   f1_per_class <- 2 * precision_per_class * recall_per_class /
#                   (precision_per_class + recall_per_class)

#   balacc <- mean(recall_per_class, na.rm = TRUE)
#   macro_f1 <- mean(f1_per_class, na.rm = TRUE)
#   macro_prec <- mean(precision_per_class, na.rm = TRUE)
#   macro_rec  <- mean(recall_per_class, na.rm = TRUE)

#   po <- acc
#   pe <- sum(rowSums(cm) / n * colSums(cm) / n)
#   kappa <- (po - pe) / (1 - pe)

# # ####################–––––––––––––––––––##########################

# cat("I: mae =", mae, "\n")
# cat("J: accuracy =", acc, "\n")
# cat("K: balanced accuracy =", balacc, "\n")
# cat("L: macro F1 =", macro_f1, "\n")
# cat("M: precision =", macro_prec, "\n")
# cat("N: recall =", macro_rec, "\n")
# cat("O: pe =", pe, "\n")
# cat("P: kappa =", kappa, "\n")

#   metrics_df <- data.frame(
#     metric = c("MAE", "Accuracy", "Balanced Accuracy", "Macro F1", "Kappa",
#                "Precision", "Recall"),
#     value  = c(mae, acc, balacc, macro_f1, kappa, macro_prec, macro_rec)
#   )


#   bin_acc  <- sum(diag(cm_bin)) / sum(cm_bin)
#   bin_rec  <- cm_bin["1", "1"] / sum(cm_bin["1", ])
#   bin_prec <- cm_bin["1", "1"] / sum(cm_bin[, "1"])
#   bin_f1   <- 2 * bin_prec * bin_rec / (bin_prec + bin_rec)
#   bin_recall_per_class <- diag(cm_bin) / rowSums(cm_bin)
#   bin_balacc <- mean(bin_recall_per_class, na.rm = TRUE)


# cat("J: bin_accuracy =", bin_acc, "\n")
# cat("K: balanced accuracy =", bin_balacc, "\n")
# cat("L: bin_ F1 =", bin_f1, "\n")
# cat("M: bin_precision =", bin_prec, "\n")
# cat("N: bin_recall =", bin_rec, "\n")

# print("Binary metrics")

#   bin_metrics_df <- data.frame(
#     metric = c("Binary Accuracy", "Binary Balanced Accuracy", "Binary F1",
#                "Binary Precision", "Binary Recall"),
#     value  = c( bin_acc, bin_balacc, bin_f1, bin_prec, bin_rec)
#   )
#   print("Ending metrics evaluation")
#   list(label = label, n_total = n,
#        metrics = metrics_df, binary = bin_metrics_df,
#        confusion = cm, confusion_bin = cm_bin,
#        dichotomize_threshold = dichotomize_threshold,
#        per_site = site_cm, sites_used = names(site_cm), sites_excluded = no_validation)
# }

# #' Convenience side-by-side comparison of two ds.semiOPBARTEvaluate() results
# #' (e.g. Architecture D vs Option C on the same held-out test rows).
# compare_semiOPBART_evaluations <- function(eval_a, eval_b) {
#   m <- merge(eval_a$metrics, eval_b$metrics, by = "metric",
#              suffixes = paste0("_", c(eval_a$label, eval_b$label)))
#   b <- merge(eval_a$binary, eval_b$binary, by = "metric",
#              suffixes = paste0("_", c(eval_a$label, eval_b$label)))
#   list(metrics = m, binary = b)
# }



# dsSemiOPBARTEvaluate.R
# ---------------------------------------------------------------------------
# Federated evaluation, mirroring evaluate_semiopbart()'s full metric set
# (MAE, accuracy, balanced accuracy, macro F1, kappa, precision, recall, plus
# the binary 0-vs-non-0 versions and both confusion matrices) -- and, as of
# this version, evaluate_semiopbart() mirrors THIS file's exact formulas
# instead of yardstick's, specifically so federated-vs-non-federated
# comparisons are apples-to-apples by construction rather than dependent on
# two different packages' internal tie-breaking rules happening to agree.
#
# KEY DESIGN POINT: every one of those metrics is a function of the
# confusion matrix alone -- and confusion-matrix cells are just counts,
# which sum safely across sites (subject to the usual nfilter floor). So
# unlike training, evaluation does NOT need an iterative round-trip
# protocol at all: one local confusion-matrix computation per site, one
# aggregation, done. This file is intentionally generic -- it works
# identically whether `newobj_pred` was produced by
# dsSemiOPBARTTrainTest.R (Architecture D) or by the predict addendum in
# dsSemiOPBARTLocalCombine.R (Option C), which is what makes the two
# directly comparable in the orchestrator.
# ---------------------------------------------------------------------------

# ---- client side ------------------------------------------------------------

#' Ordinal mean absolute error computed DIRECTLY from an aggregated
#' (truth x predicted) confusion matrix, i.e. from counts that are already
#' safe to have crossed the DataSHIELD boundary (see ds.semiOPBARTEvaluate()
#' / semiOPBARTLocalConfusionDS()) -- no new disclosure surface, and no
#' reliance on raw per-patient predictions ever reaching the client (they
#' never do, and never should).
#'
#' cm's rows and columns are both ordered by `levels` (see
#' semiOPBARTLocalConfusionDS(): `factor(..., levels = levels)` for both
#' truth and prediction), so `outer(lv, lv, "-")` lines up with cm exactly.
#' `levels` here must match cm's actual dimensions -- pass the MULTICLASS
#' level set for the multiclass cm, and c("0","1") for a binary cm, not the
#' same `levels` for both (see .semiOPBART_ordinalMAE()'s two call sites in
#' ds.semiOPBARTEvaluate() below for why this distinction matters).
.semiOPBART_ordinalMAE <- function(cm, levels) {
  lv <- suppressWarnings(as.numeric(as.character(levels)))
  if (anyNA(lv)) lv <- seq_along(levels) - 1  # fallback: treat as equally-spaced ranks
  dist_mat <- abs(outer(lv, lv, "-"))
  sum(cm * dist_mat) / sum(cm)
}

#' Aggregate confusion matrices across sites and compute the full
#' evaluate_semiopbart()-equivalent metric set from the pooled matrix.
#'
#' @param pred_obj  name of the local prediction object at each site
#' @param levels    full ordinal level set (same as used in prepare/predict)
#' @param label     a label for this evaluation run (e.g. "Architecture D",
#'                  "Option C"), carried through for the comparison step
#' @param validation_col, validation_data.name  see
#'   semiOPBARTLocalConfusionDS(). If validation_col is given, sites
#'   without that column (or without enough non-missing values) are simply
#'   excluded from the aggregate, with a message naming them, rather than
#'   failing the whole evaluation.
#' @param dichotomize_threshold  NULL (default): binary evaluation uses the
#'   original baseline-vs-rest split. A number (e.g. 2): binary evaluation
#'   dichotomizes at this threshold instead (truth/prediction >= threshold
#'   -> "1"). See semiOPBARTLocalConfusionDS() for the full rationale.
#' @export
ds.semiOPBARTEvaluate <- function(pred_obj = "semiOPBART_pred", levels,
                                   nfilter = 5, datasources = NULL,
                                   label = "model",
                                   validation_col = NULL,
                                   validation_data.name = NULL,
                                   dichotomize_threshold = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()

  site_cm_raw <- DSI::datashield.aggregate(datasources,
      call("semiOPBARTLocalConfusionDS", pred_obj, semiOPBART_toSerialize(levels), nfilter,
           if (is.null(validation_col)) "null" else validation_col,
           if (is.null(validation_data.name)) "null" else validation_data.name,
           if (is.null(dichotomize_threshold)) "null" else as.character(dichotomize_threshold)))

  no_validation <- names(site_cm_raw)[vapply(site_cm_raw, is.null, logical(1))]
  if (length(no_validation) && !is.null(validation_col))
    message("ds.semiOPBARTEvaluate(): no '", validation_col, "' validation ",
            "available at: ", paste(no_validation, collapse = ", "),
            " -- excluded from this evaluation, not treated as an error.")

  # FIX: an earlier version of this file re-assigned `site_cm <-
  # site_cm_raw` right after this filtering step, silently discarding it
  # (every NULL site would then hit `NULL[["confusion"]]` -> NULL inside
  # the Reduce() below, propagating to numeric(0) and quietly wrecking the
  # whole aggregation with no error). Filtered result kept, not overwritten.
  site_cm <- site_cm_raw[!vapply(site_cm_raw, is.null, logical(1))]
  if (length(site_cm) == 0)
    stop("no site could be validated -- ",
         if (!is.null(validation_col))
           paste0("'", validation_col, "' was unavailable everywhere")
         else "no predictions with a training-label truth were found")

  # Aggregation uses the UNMASKED $confusion (not $confusion_reported) so
  # cross-site sums are exact even when an individual site's cell was too
  # small to report on its own -- the pooled cell is what actually needs
  # the nfilter check, not each site's contribution to it.
  cm     <- Reduce(`+`, lapply(site_cm, `[[`, "confusion"))
  cm_bin <- Reduce(`+`, lapply(site_cm, `[[`, "confusion_bin"))
  if (any(cm > 0 & cm < nfilter))
    warning("a pooled multiclass confusion cell is below nfilter=", nfilter,
            " even after aggregation -- treat that cell's contribution to ",
            "per-class metrics with caution")

  n <- sum(cm)
  acc <- sum(diag(cm)) / n
  mae <- .semiOPBART_ordinalMAE(cm, levels)

  recall_per_class <- diag(cm) / rowSums(cm)
  precision_per_class <- diag(cm) / colSums(cm)
  f1_per_class <- 2 * precision_per_class * recall_per_class /
                  (precision_per_class + recall_per_class)

  balacc <- mean(recall_per_class, na.rm = TRUE)
  macro_f1 <- mean(f1_per_class, na.rm = TRUE)
  macro_prec <- mean(precision_per_class, na.rm = TRUE)
  macro_rec  <- mean(recall_per_class, na.rm = TRUE)

  po <- acc
  pe <- sum(rowSums(cm) / n * colSums(cm) / n)
  kappa <- (po - pe) / (1 - pe)
cat("mae =", mae, "\n")
cat("accuracy =", acc, "\n")
cat("balanced accuracy =", balacc, "\n")
cat("macro F1 =", macro_f1, "\n")
cat("precision =", macro_prec, "\n")
cat("recall =", macro_rec, "\n")
cat("pe =", pe, "\n")
cat("P: kappa =", kappa, "\n")
  metrics_df <- data.frame(
    metric = c("MAE", "Accuracy", "Balanced Accuracy", "Macro F1", "Kappa",
               "Precision", "Recall"),
    value  = c(mae, acc, balacc, macro_f1, kappa, macro_prec, macro_rec)
  )

  # FIX: bin_mae previously passed the FULL multiclass `levels` (e.g.
  # length 5) to .semiOPBART_ordinalMAE() alongside a 2x2 BINARY cm --
  # outer(lv, lv, "-") with a length-5 lv gives a 5x5 dist_mat, and
  # `cm_bin * dist_mat` (2x2 * 5x5) is non-conformable, a hard error in R.
  # Binary confusion matrices are always exactly c("0","1") (see
  # semiOPBARTLocalConfusionDS()), so that's what gets passed here, not
  # the multiclass `levels` argument.
  bin_mae  <- .semiOPBART_ordinalMAE(cm_bin, c("0", "1"))
  bin_acc  <- sum(diag(cm_bin)) / sum(cm_bin)
  bin_rec  <- cm_bin["1", "1"] / sum(cm_bin["1", ])
  bin_prec <- cm_bin["1", "1"] / sum(cm_bin[, "1"])
  bin_f1   <- 2 * bin_prec * bin_rec / (bin_prec + bin_rec)
  bin_recall_per_class <- diag(cm_bin) / rowSums(cm_bin)
  bin_balacc <- mean(bin_recall_per_class, na.rm = TRUE)

cat("bin_mae =", bin_mae, "\n")
cat("bin_accuracy =", bin_acc, "\n")
cat("bin_balanced accuracy =", bin_balacc, "\n")
cat("bin_ F1 =", bin_f1, "\n")
cat("bin_precision =", bin_prec, "\n")
cat("bin_recall =", bin_rec, "\n")
  bin_metrics_df <- data.frame(
    metric = c("Binary MAE", "Binary Accuracy", "Binary Balanced Accuracy", "Binary F1",
               "Binary Precision", "Binary Recall"),
    value  = c(bin_mae, bin_acc, bin_balacc, bin_f1, bin_prec, bin_rec)
  )

  list(label = label, n_total = n,
       metrics = metrics_df, binary = bin_metrics_df,
       confusion = cm, confusion_bin = cm_bin,
       dichotomize_threshold = dichotomize_threshold,
       per_site = site_cm, sites_used = names(site_cm), sites_excluded = no_validation)
}

#' Convenience side-by-side comparison of two ds.semiOPBARTEvaluate() results
#' (e.g. Architecture D vs Option C on the same held-out test rows).
#' @export
compare_semiOPBART_evaluations <- function(eval_a, eval_b) {
  m <- merge(eval_a$metrics, eval_b$metrics, by = "metric",
             suffixes = paste0("_", c(eval_a$label, eval_b$label)))
  b <- merge(eval_a$binary, eval_b$binary, by = "metric",
             suffixes = paste0("_", c(eval_a$label, eval_b$label)))
  list(metrics = m, binary = b)
}
