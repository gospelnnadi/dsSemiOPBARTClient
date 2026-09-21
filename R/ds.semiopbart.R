# dsSemiOPBARTClient.R
# ---------------------------------------------------------------------------
# Client-side orchestrator for federated semi-OPBART, Architecture D.
# Every outbound call() argument beyond a bare scalar is Serialize-encoded here
# via semiOPBART_toSerialize() -- see dsSemiOPBARTSerialize.R for why, and
# dsSemiOPBARTBase.R for the matching server-side decode calls.
# ---------------------------------------------------------------------------

#source("ds.serialize.R")

ds.semiOPBARTTrainE <- function(data.name = "D_train", datasources = NULL,
                                 num_tree = 20, k = 1,
                                 num_burn = 1000, num_save = 1000,
                                 sync_every = 5,
                                 warmup_burn = 100,
                                 nfilter_threshold = 5,
                                 threshold_method =  c("exact","percentile"),
                                 state_name = ".semiOPBART_state", seed = 35,
                                 quiet = TRUE,
                                   diagnose = TRUE,
                                   diagnose_every = 50,
                                   theta_alarm = 20,
                                   us_alarm = 50,
                                    sd = 1) {
  set.seed(seed)
  threshold_method <- match.arg(threshold_method)
  
  # ---- FILTER: Only connections with train object ----
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  
  # Use the shared filter
  trainable_conns <- .filter_trainable_connections(
    conns = datasources,
    train.name = data.name,
    verbose = TRUE
  )
  print(paste0("[ds.semiOPBARTTrainE] Trainable connections: ", paste(names(trainable_conns), collapse = ", ")))
  if (length(trainable_conns) == 0) {
    stop("No trainable sites. All sites were auto-demoted.")
  }
  
  site_names <- names(trainable_conns)
  
  # PERFORMANCE: suppress progress bars if quiet
  if (quiet) {
    old_progress <- getOption("datashield.progress")
    options(datashield.progress = FALSE)
    on.exit(options(datashield.progress = old_progress), add = TRUE)
  }

  site_meta <- DSI::datashield.aggregate(trainable_conns,
      call("semiOPBARTLocalInitEDS", data.name, num_tree, k, nfilter_threshold, state_name, seed = seed))
  
  p <- site_meta[[1]]$p
  J <- site_meta[[1]]$J
  if (!all(vapply(site_meta, function(m) m$p == p && m$J == J, logical(1))))
    stop("sites disagree on p (ncol W) or J (number of categories)")


# ADD THIS:
site_n <- setNames(vapply(site_meta, function(m) as.integer(m$n), integer(1)),
                    names(trainable_conns))

  theta <- rep(1, p)
  us <- 0:(J - 2); us[1] <- 0
  message(sprintf("[ds.semiOPBARTTrainE] Training on %d site(s): %s", 
                  length(site_names), paste(site_names, collapse = ", ")))
  message(sprintf("[ds.semiOPBARTTrainE] Initialized with p=%d, J=%d", p, J))

  N <- num_burn + num_save
  theta_draws <- matrix(NA_real_, N, p)
  us_draws <- matrix(NA_real_, N, length(us))

  total_sweeps_done <- 0

  while (total_sweeps_done < N) {

    sweeps_remaining <- N - total_sweeps_done
    num_local_sweeps <- min(sync_every, sweeps_remaining)
    is_last_block <- (total_sweeps_done + num_local_sweeps >= N)

    # -----------------------------------------------------------------
    # 1. Run local MCMC on EVERY site - NO client communication
    #    (Architecture E: no tree swapping ever)
    # -----------------------------------------------------------------
    message(sprintf("[ds.semiOPBARTTrainE] Running %d local sweeps (total %d/%d)", num_local_sweeps, total_sweeps_done + num_local_sweeps, N))
    site_results <- DSI::datashield.aggregate(datasources,
        call("semiOPBARTLocalMCMC",
             num_local_sweeps,
             semiOPBART_toSerialize(theta),
             semiOPBART_toSerialize(us),
             FALSE,  # No tree swapping
             0,
             "null",
             "null",
             nfilter_threshold,
             J,
             threshold_method,
             state_name,
             seed = seed,
              sd = sd))
    message(sprintf("[ds.semiOPBARTTrainE] Completed %d local sweeps (total %d/%d)", num_local_sweeps, total_sweeps_done + num_local_sweeps, N))
    total_sweeps_done <- total_sweeps_done + num_local_sweeps

    # -----------------------------------------------------------------
    # 2. AGGREGATE theta across sites
    # -----------------------------------------------------------------
    WtW <- Reduce(`+`, lapply(site_results, `[[`, "WtW"))
    WtZr <- Reduce(`+`, lapply(site_results, `[[`, "WtZr"))

    theta_hat <- solve(WtW, WtZr)
    theta_sigma <- solve(WtW)
    theta <- as.numeric(mvtnorm::rmvnorm(1, theta_hat, theta_sigma))

    # -----------------------------------------------------------------
    # 3. AGGREGATE thresholds across sites
    # -----------------------------------------------------------------
    th_stats <- lapply(site_results, `[[`, "th_stats")
     print(paste0 ("LocalMCMC: length of theta_stats ", length(theta_hat) ))

    us <- update_thresholds_from_site_stats(th_stats, us)
    message(sprintf("[ds.semiOPBARTTrainE] After aggregation: max|theta|=%.3g, max|us|=%.3g", max(abs(theta)), max(abs(us))))
      # DIAGNOSTIC (read-only -- does not change theta/us/WtW at all, or
    # whether/when a downstream do_gibbs() chol()/inv_sympd() failure
    # happens): periodically, and immediately on either theta or us
    # crossing an "this looks like divergence, not real data" magnitude,
    # report where things stand. The point is specifically to tell apart
    # "unstable from sweep 1" (near-collinear w_features -- a data
    # problem) from "fine for a while, then drifts" (unbounded
    # identifiability under the flat prior -- see this function's
    # docstring) the NEXT time this fails, instead of only ever seeing
    # the crash itself with no trajectory leading up to it.
    if (diagnose) {
      theta_extreme <- max(abs(theta)) > theta_alarm
      us_extreme <- max(abs(us)) > us_alarm
      if (theta_extreme || us_extreme || total_sweeps_done %% diagnose_every == 0 || total_sweeps_done == 1) {
        kappa_WtW <- tryCatch(kappa(WtW, exact = FALSE), error = function(e) NA_real_)
        message(sprintf(
          "[ds.semiOPBARTTrainE] sweep %d/%d: max|theta|=%.3g, max|us|=%.3g, kappa(WtW)=%.3g%s",
          total_sweeps_done, N, max(abs(theta)), max(abs(us)), kappa_WtW,
          if (theta_extreme || us_extreme) "  <-- EXTREME, likely divergence starting here" else ""))
      }
    }

    # Store draws for this block
    block_start <- total_sweeps_done - num_local_sweeps + 1
    block_end <- total_sweeps_done

    for (idx in block_start:block_end) {
      theta_draws[idx, ] <- theta
      us_draws[idx, ] <- us
    }
    
  }

  keep <- (num_burn + 1):N
  list(theta_draws = theta_draws[keep, , drop = FALSE],
       us_draws = us_draws[keep, , drop = FALSE],
       theta_mean = colMeans(theta_draws[keep, , drop = FALSE]),
       us_mean = colMeans(us_draws[keep, , drop = FALSE]),
       site_names = site_names, data.name = data.name,
       num_tree = num_tree, k = k, site_n = site_n, state_name = state_name)
}






# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


#' Combine per-site threshold statistics into the next global threshold
#' vector. CONFIRMED against the real smopbart() source
#' (`us_p[i,j] = runif(1, min = max(Z[Y==j], us_p[i,j-1]), max =
#' min(Z[Y==j+1], us_p[i-1,j+1]))`, or Inf instead of that last term for
#' the final threshold): threshold j's LOWER bound comes from category
#' j's UPPER tail (this function's `hi` for category j), and its UPPER
#' bound comes from category (j+1)'s LOWER tail (this function's `lo` for
#' category j+1) -- NOT the same category's own lo/hi for both bounds,
#' which is what an earlier version of this function did. That earlier
#' version was a real bug, not just a disclosure-motivated approximation:
#' it would have pulled thresholds toward the wrong values regardless of
#' the percentile-vs-exact-min/max disclosure substitution, which is a
#' separate, intentional approximation (percentiles standing in for the
#' real algorithm's exact max(Z[Y==j])/min(Z[Y==j+1])).
#'
#' THE LAST THRESHOLD IS NOT A SPECIAL CASE: its upper bound is the TOP
#' category's own `lo` statistic (category J's 5th percentile of Z) --
#' the EXACT SAME quantity every interior threshold already uses from its
#' own upper-neighbor category, since semiOPBARTLocalThresholdStatsDS()
#' already reports all J categories, including the last one. An earlier
#' version of this function used upper_bound = Inf for the last
#' threshold (no real bound at all), and a since-reverted local patch
#' replaced that with `lower_bound + abs(rnorm(1))` -- an always-positive
#' perturbation with nothing pulling it back down, which random-walks
#' upward without limit sweep after sweep (visible as the last threshold
#' climbing into the hundreds after enough iterations). Neither was a
#' real posterior draw; both are removed here in favor of the same
#' bounded-uniform draw every other threshold already uses.
#'
#' @param th_stats  per-site output of semiOPBARTLocalThresholdStatsDS(),
#'   one entry per category "1".."J" (not per threshold)
#' @param us_prev   previous sweep's threshold vector, length J-1, us[1]
#'   always fixed at 0 (identifiability restriction, unchanged here)
update_thresholds_from_site_stats <- function(th_stats, us_prev) {
  J <- length(th_stats[[1]])           # categories, not thresholds
  n_thresh <- length(us_prev)          # == J - 1
  new_us <- us_prev

  get_stat <- function(cat, which)
    vapply(th_stats, function(s) s[[as.character(cat)]][[which]], numeric(1))

  for (j in 2:n_thresh) {              # matches real smopbart()'s
                                        # `for(j in 2:ncol(us_p))`; column
                                        # 1 stays fixed at 0
    hi_j  <- get_stat(j,     "hi")     # this threshold's LOWER category
    lo_j1 <- get_stat(j + 1, "lo")     # this threshold's UPPER category --
                                        # ALWAYS finite and available, even
                                        # when j+1 == J (the top category):
                                        # see the note above on why there is
                                        # no "last threshold" special case

    # sequential (systematic-scan) Gibbs: use new_us[j-1], the value THIS
    # SAME sweep just drew, not the previous sweep's stale us_prev[j-1] --
    # keeps every draw consistent with what was just sampled immediately
    # before it, rather than only self-consistent one sweep later.
    lower_bound <- max(c(hi_j[!is.na(hi_j)], new_us[j - 1]))
    upper_bound <- min(c(lo_j1[!is.na(lo_j1)],
                         if (j < n_thresh) us_prev[j + 1] else Inf))

    # guard against a degenerate/inverted interval (possible with sparse
    # categories after nfilter suppression at some sites, or an unlucky
    # sweep) -- fall back to holding the previous value rather than
    # erroring or drawing from an inverted/empty range. Inf can only
    # still appear here if EVERY site's category-J "lo" was NA (i.e. that
    # whole category was too sparse everywhere this sweep) -- a real,
    # rare data condition, not a bug; holding us_prev[j] is the right
    # fallback for it too.
    new_us[j] <- if (is.finite(lower_bound) && is.finite(upper_bound) &&
                      upper_bound > lower_bound)
      stats::runif(1, lower_bound, upper_bound) else us_prev[j]
  }
  new_us[1] <- 0
  new_us
}



#' Variable importance for a trained Architecture D or E fit: how often
#' each X feature is used as a splitting variable, pooled (n-weighted, by
#' default) across every trained site's CURRENT forest via the new
#' semiOPBARTLocalVarCountsDS() aggregate (semiopbartDS.R) -- a snapshot
#' of the final trained forest, not a full posterior (see that function's
#' docstring for why). For Architecture F, use
#' semiOPBART_explainVarImportanceF() instead (dsSemiOPBARTExplain.R) --
#' F's local smopbart() already exposes this directly, no extra
#' aggregate() round trip needed.
#'
#' @param fit  a ds.semiOPBARTTrain()/ds.semiOPBARTTrainE() result --
#'   needs $site_names and $state_name, both of which those return
#' @param datasources  defaults to every connection in fit$site_names;
#'   pass an explicit subset to restrict to sites you know are still
#'   trained/connected
#' @param weight_by_n  n-weighted average across sites (default) or
#'   equal-weighted
#' @export
ds.semiOPBARTVarImportance <- function(fit, datasources = NULL, weight_by_n = TRUE) {
  if (is.null(fit$site_names) || is.null(fit$state_name))
    stop("ds.semiOPBARTVarImportance(): `fit` needs both $site_names and ",
         "$state_name -- pass the object returned by ds.semiOPBARTTrain() ",
         "or ds.semiOPBARTTrainE(), not a predict/evaluate result.")
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  conns <- datasources[fit$site_names]

  site_results <- DSI::datashield.aggregate(conns,
      call("semiOPBARTLocalVarCountsDS", fit$state_name))

  n <- vapply(site_results, `[[`, numeric(1), "n")
  w <- if (weight_by_n) n / sum(n) else rep(1 / length(site_results), length(site_results))

  counts_list <- lapply(site_results, `[[`, "var_counts")
  lens <- vapply(counts_list, length, integer(1))
  if (length(unique(lens)) > 1)
    stop("ds.semiOPBARTVarImportance(): sites disagree on the number of X ",
         "features (", paste(unique(lens), collapse = ", "), ") -- did ",
         "they all train on the same formula?")

  pooled <- Reduce(`+`, Map(function(cnt, wi) wi * cnt, counts_list, w))
  names(pooled) <- names(counts_list[[1]])   # NULL if the server-side X had no colnames
  pooled
}



# ds.semiopbart.R - Add this helper function

#' Filter connections to only those with a trainable object
#' 
#' This is the single shared choke point for all trainers:
#' - Checks each connection for the presence of `train.name`
#' - Returns only connections where the object exists
#' - Logs which sites are excluded (with reason)
#' - No data is ever returned to the client
#' 
#' @param conns  list of DataSHIELD connections
#' @param train.name  name of the training object to check for
#' @param verbose  print exclusion messages
#' @return list of connections that have the train object
#' @export
.filter_trainable_connections <- function(conns, train.name = "D_train", 
                                           verbose = TRUE) {
  if (is.null(conns)) {
    conns <- DSI::datashield.connections_find()
  }
  
  site_names <- names(conns)
  if (is.null(site_names)) {
    stop(".filter_trainable_connections: datasources must be named")
  }
  
  # Check each site for the train object
  has_train <- logical(length(conns))
  names(has_train) <- site_names
  
  for (i in seq_along(conns)) {
    result <- tryCatch({
      DSI::datashield.aggregate(
        conns[i],
        call("semiOPBARTLocalCheckExistsDS", train.name)
      )[[1]]
    }, error = function(e) {
      # If the server function itself fails, treat as FALSE
      FALSE
    })
    has_train[i] <- isTRUE(result)
  }
  
  # Filter to only sites with train object
  trainable_idx <- which(has_train)
  trainable_conns <- conns[trainable_idx]
  
  # Report excluded sites
  excluded <- site_names[!has_train]
  if (length(excluded) > 0 && verbose) {
    message(sprintf(
      ".filter_trainable_connections: Excluded %d site(s) with no '%s': %s",
      length(excluded), train.name, paste(excluded, collapse = ", ")
    ))
    message("  (These sites were auto-demoted by ds.semiOPBARTSplit() - they have too few labeled rows to train)")
  }
  
  if (length(trainable_conns) == 0) {
    stop(sprintf(
      ".filter_trainable_connections: No sites have '%s'. All sites were auto-demoted.",
      train.name
    ))
  }
  
  trainable_conns
}