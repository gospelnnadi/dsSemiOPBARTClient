# dsSemiOPBARTClient.R
# ---------------------------------------------------------------------------
# Client-side orchestrator for federated semi-OPBART, Architecture D.
# Every outbound call() argument beyond a bare scalar is Serialize-encoded here
# via semiOPBART_toSerialize() -- see dsSemiOPBARTSerialize.R for why, and
# dsSemiOPBARTBase.R for the matching server-side decode calls.
# ---------------------------------------------------------------------------

#source("ds.serialize.R")
# dsSemiOPBARTClient.R - Updated Architecture D with batched communication

# dsSemiOPBARTClient.R - Updated with batched communication

#' @param data.name    the TRAIN object from ds.semiOPBARTSplit()
#' @param sync_every  perform federated aggregation and tree exchange every
#'   `sync_every` sweeps. Between these, sites run local MCMC ENTIRELY
#'   on the server with NO client communication.
#' @param n_swap      how many trees each site exchanges with its swap partner
#'   each swap event. Must be <= num_tree.
#' @param threshold_method  "percentile" (default, disclosure-safe) or "exact"
#' @export
ds.semiOPBART <- function(data.name = "semiOPBART_train", datasources = NULL,
                           num_tree = 20, k = 1,
                           num_burn = 1000, num_save = 1000,
                           sync_every = 50, n_swap = 4,
                           nfilter_threshold = 5,
                           threshold_method = c("percentile", "exact"),
                           state_name = ".semiOPBART_state",
                           allow_local_ecdf = FALSE, seed = 35) {
  set.seed(seed)
  threshold_method <- match.arg(threshold_method)
  if (n_swap > num_tree)
    stop("n_swap (", n_swap, ") cannot exceed num_tree (", num_tree, ")")

  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  site_names <- names(datasources)
  n_sites <- length(site_names)

  # Initialize EVERY site's local forest
  site_meta <- DSI::datashield.aggregate(datasources,
      call("semiOPBARTLocalInitDS", data.name, num_tree, k, 5, state_name, allow_local_ecdf, seed = seed))

  p <- site_meta[[1]]$p
  J <- site_meta[[1]]$J
  if (!all(vapply(site_meta, function(m) m$p == p && m$J == J, logical(1))))
    stop("sites disagree on p (ncol W) or J (number of categories)")

  theta <- rep(0, p)
  us <- 0:(J - 2); us[1] <- 0

  N <- num_burn + num_save
  theta_draws <- matrix(NA_real_, N, p)
  us_draws <- matrix(NA_real_, N, length(us))
  
  # Track where we are in the MCMC
  total_sweeps_done <- 0
  
  # We'll run in blocks of sync_every local sweeps
  while (total_sweeps_done < N) {
    
    # Determine how many local sweeps to run before next communication
    sweeps_remaining <- N - total_sweeps_done
    num_local_sweeps <- min(sync_every, sweeps_remaining)
    is_last_block <- (total_sweeps_done + num_local_sweeps >= N)
    
    # Check if this block should include tree swapping
    # Swap happens at the END of a block (after the local sweeps)
    perform_swap <- (n_sites > 1 && 
                     num_local_sweeps == sync_every && 
                     !is_last_block)
    
    # -----------------------------------------------------------------
    # 1. Run local MCMC on EVERY site - NO client communication in between!
    # -----------------------------------------------------------------
    site_results <- DSI::datashield.aggregate(datasources,
        call("semiOPBARTLocalMCMC", 
             num_local_sweeps,
             semiOPBART_toSerialize(theta),
             semiOPBART_toSerialize(us),
             perform_swap,  # This only tells the server to prepare for swap
             n_swap,
             "null",  # swap_targets will be set below after we know indices
             "null",  # swap_incoming will be set below
             nfilter_threshold,
             J,
             threshold_method,
             state_name,
             seed = seed))
    
    # Update total sweeps done
    total_sweeps_done <- total_sweeps_done + num_local_sweeps
    
    # -----------------------------------------------------------------
    # 2. AGGREGATE theta across sites
    # -----------------------------------------------------------------
    WtW <- Reduce(`+`, lapply(site_results, `[[`, "WtW"))
    WtZr <- Reduce(`+`, lapply(site_results, `[[`, "WtZr"))
    theta_hat <- solve(WtW) %*% WtZr
    theta_sigma <- solve(WtW)
    theta <- as.numeric(mvtnorm::rmvnorm(1, theta_hat, theta_sigma))
    
    # -----------------------------------------------------------------
    # 3. AGGREGATE thresholds across sites
    # -----------------------------------------------------------------
    th_stats <- lapply(site_results, `[[`, "th_stats")
    us <- update_thresholds_from_site_stats(th_stats, us)
    
    # Store draws for this block
    block_start <- total_sweeps_done - num_local_sweeps + 1
    block_end <- total_sweeps_done
    
    # We use the NEW theta/us for ALL draws in this block
    # (This is a design choice - using the post-sync values)
    for (idx in block_start:block_end) {
      theta_draws[idx, ] <- theta
      us_draws[idx, ] <- us
    }
    
    # -----------------------------------------------------------------
    # 4. TREE SWAPPING (if applicable)
    # -----------------------------------------------------------------
    if (perform_swap && n_sites > 1) {
      
      # Choose which trees to swap at each site
      export_idx <- setNames(
        lapply(site_names, function(s) sort(sample(0:(num_tree - 1), n_swap))),
        site_names)
      
      # Decode trees from each site
      export_trees <- setNames(lapply(site_names, function(s) {
        all_trees <- semiOPBART_fromSerialize(site_results[[s]]$trees_Serialize)
        all_trees[export_idx[[s]] + 1]
      }), site_names)
      
      # Perform swaps
      # FIX: was `for (k in seq_len(n_sites))`, silently shadowing this
      # function's own `k` parameter (BART shrinkage) -- see the identical
      # fix in ds.semiOPBARTWarmupBlock().
      for (swap_i in seq_len(n_sites)) {
        sender <- site_names[swap_i]
        receiver <- site_names[if (swap_i == n_sites) 1 else swap_i + 1]
        
        DSI::datashield.aggregate(datasources[receiver],
            call("semiOPBARTLocalSwapDS", 
                 semiOPBART_toSerialize(export_idx[[receiver]]),
                 semiOPBART_toSerialize(export_trees[[sender]]),
                 state_name,
                 seed = seed))
      }
    }
    
    # After swap, broadcast updated theta/us to ALL sites for next block
    # (The swap function doesn't change theta/us, but we need the new values
    #  for the next block's local sweeps)
    # Actually, theta/us are already in the client environment.
    # For the next block, we'll pass the updated theta/us to LocalMCMC.
    
    # But wait - the sites still have the OLD theta/us in their local state!
    # We need to update them before the next block.
    # We'll do this implicitly by passing the new theta/us to the next 
    # LocalMCMC call. The sites will use these as starting values.
    
    # For the draws we just stored, we used the NEW theta/us.
    # For the next block, we'll start from these new values.
  }

  keep <- (num_burn + 1):N
  list(theta = theta_draws[keep, , drop = FALSE],
       us = us_draws[keep, , drop = FALSE],
       data.name = data.name,
       num_tree = num_tree, k = k,
       site_names = site_names, state_name = state_name)
}


# #' @export
# ds.semiOPBARTTrain <- ds.semiOPBART   # explicit train-only alias, see
#                                        # dsSemiOPBARTTrainTest.R


                                       
                                                                          
# dsSemiOPBARTClient.R - Architecture E with batched communication

#' @param data.name  the TRAIN object at each site
#' @param sync_every  perform federated theta/us aggregation every
#'   `sync_every` sweeps. Between these, sites run local MCMC ENTIRELY
#'   on the server with NO client communication.
#' @export
ds.semiOPBARTTrainE_old <- function(data.name = "D_train", datasources = NULL,
                                 num_tree = 20, k = 1,
                                 num_burn = 1000, num_save = 1000,
                                 sync_every = 50,
                                 nfilter_threshold = 5,
                                 threshold_method = c("percentile", "exact"),
                                 state_name = ".semiOPBART_state", seed = 35) {
  set.seed(seed)                                
  threshold_method <- match.arg(threshold_method)
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  site_names <- names(datasources)

  site_meta <- DSI::datashield.aggregate(datasources,
      call("semiOPBARTLocalInitEDS", data.name, num_tree, k, 5, state_name, seed = seed))
  p <- site_meta[[1]]$p
  J <- site_meta[[1]]$J
  if (!all(vapply(site_meta, function(m) m$p == p && m$J == J, logical(1))))
    stop("sites disagree on p (ncol W) or J (number of categories)")

  theta <- rep(0, p)
  us <- 0:(J - 2); us[1] <- 0

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
             seed = seed))
    
    total_sweeps_done <- total_sweeps_done + num_local_sweeps
    
    # -----------------------------------------------------------------
    # 2. AGGREGATE theta across sites
    # -----------------------------------------------------------------
    WtW <- Reduce(`+`, lapply(site_results, `[[`, "WtW"))
    WtZr <- Reduce(`+`, lapply(site_results, `[[`, "WtZr"))
    theta_hat <- solve(WtW) %*% WtZr
    theta_sigma <- solve(WtW)
    theta <- as.numeric(mvtnorm::rmvnorm(1, theta_hat, theta_sigma))
    
    # -----------------------------------------------------------------
    # 3. AGGREGATE thresholds across sites
    # -----------------------------------------------------------------
    th_stats <- lapply(site_results, `[[`, "th_stats")
    us <- update_thresholds_from_site_stats(th_stats, us)
    
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
       num_tree = num_tree, k = k, state_name = state_name)
}
# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Ring-topology tree swap: site k's chosen n_swap trees go to site k+1
#' (mod n_sites), landing in site (k+1)'s OWN chosen n_swap slots -- i.e.
#' every site simultaneously gives away n_swap trees and receives n_swap
#' different ones INTO THE SAME local slot positions it gave up, so each
#' site's total stays exactly num_tree, always (a strict swap, never a
#' net transfer). Both sides' chosen indices are picked CLIENT-SIDE (pure
#' index bookkeeping, no data involved, no server round-trip needed just
#' to choose them).
#'
#' @param grow_out   this sweep's semiOPBARTLocalGrowDS() output, keyed by
#'   site -- each site's trees_Serialize covers its FULL current num_tree-tree
#'   forest (not a subset), so any local index can be sliced out of it
#' @param num_tree   this (and every) site's own forest size
#' @param n_swap     how many trees to exchange per site per swap event
swap_trees <- function(datasources, site_names, grow_out, num_tree, n_swap,
                       state_name = ".semiOPBART_state",seed = 35) {
  set.seed(seed)
  n_sites <- length(site_names)
  export_idx <- setNames(
    lapply(site_names, function(s) sort(sample(0:(num_tree - 1), n_swap))),
    site_names)

  # decode each site's full just-grown forest ONCE, slice out the trees
  # at THAT site's own export_idx (R lists are 1-indexed, tree slots are
  # 0-indexed, hence the +1)
  export_trees <- setNames(lapply(site_names, function(s) {
    all_trees <- semiOPBART_fromSerialize(grow_out[[s]]$trees_Serialize)
    all_trees[export_idx[[s]] + 1]
  }), site_names)

  for (k in seq_len(n_sites)) {
    sender   <- site_names[k]
    receiver <- site_names[if (k == n_sites) 1 else k + 1]
    incoming_Serialize   <- semiOPBART_toSerialize(unname(export_trees[[sender]]))
    target_idx_Serialize <- semiOPBART_toSerialize(export_idx[[receiver]])
    DSI::datashield.aggregate(datasources[receiver],
        call("semiOPBARTLocalSwapDS", target_idx_Serialize, incoming_Serialize, state_name, seed = seed))
  }
  invisible(NULL)
}

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



# dsSemiOPBARTClient.R - Updated with all improvements

#' Adaptive communication scheduler
#' 
#' @param theta_history  list of theta changes from each site
#' @param current_every  current communication frequency
#' @param min_every      minimum frequency (most communication)
#' @param max_every      maximum frequency (least communication)
#' @param target_change  target relative change threshold
#' @param momentum       smoothing parameter (0-1)
adaptive_scheduler <- function(theta_history, 
                               current_every,
                               min_every = 5,
                               max_every = 200,
                               target_change = 0.01,
                               momentum = 0.3) {
  
  # Average change across sites
  avg_change <- mean(sapply(theta_history, `[[`, "theta_change"), na.rm = TRUE)
  
  # If no change data, maintain current
  if (is.na(avg_change) || avg_change == 0) {
    return(current_every)
  }
  
  # Adjust frequency based on change
  # Higher change -> lower every (more frequent communication)
  ratio <- target_change / avg_change
  new_every <- current_every * ratio
  
  # Apply bounds
  new_every <- max(min_every, min(max_every, new_every))
  
  # Apply momentum for smooth transitions
  new_every <- momentum * new_every + (1 - momentum) * current_every
  
  return(round(new_every))
}


####################################################improvements#######################
#' ds.semiOPBARTImproved - Architecture D with optional adaptive communication
#' 
#' @param data.name      TRAIN object name
#' @param datasources    list of DataSHIELD connections
#' @param num_tree       number of trees per site
#' @param k              BART shrinkage parameter
#' @param num_burn       number of burn-in iterations
#' @param num_save       number of saving iterations
#' @param theta_every    frequency for theta aggregation. If NULL (default), 
#'                       uses adaptive communication. If numeric, uses fixed 
#'                       periodic communication.
#' @param us_every       frequency for threshold aggregation. If NULL (default), 
#'                       uses adaptive communication. If numeric, uses fixed 
#'                       periodic communication.
#' @param sync_every     frequency for tree swapping. If NULL (default), 
#'                       uses adaptive communication. If numeric, uses fixed 
#'                       periodic communication.
#' @param initial_every  initial communication frequency (only used when 
#'                       adaptive is enabled for a parameter)
#' @param min_every      minimum communication interval (adaptive only)
#' @param max_every      maximum communication interval (adaptive only)
#' @param warmup_burn    number of warmup iterations with aggressive communication
#' @param adapt_rate     local adaptation rate during long gaps (0-1)
#' @param n_swap         number of trees to swap per event
#' @param threshold_method "percentile" or "exact"
#' @param allow_local_ecdf allow local_ecdf with tree sharing
#' @export
ds.semiOPBARTImproved <- function(data.name = "semiOPBART_train", 
                                   datasources = NULL,
                                   num_tree = 20, 
                                   k = 1,
                                   num_burn = 1000, 
                                   num_save = 1000,
                                   theta_every = NULL,
                                   us_every = NULL,
                                   sync_every = NULL,
                                   initial_every = 50,
                                   min_every = 5,
                                   max_every = 200,
                                   warmup_burn = 100,
                                   adapt_rate = 0.1,
                                   n_swap = 4,
                                   nfilter_threshold = 5,
                                   threshold_method = c("percentile", "exact"),
                                   state_name = ".semiOPBART_state",
                                   allow_local_ecdf = TRUE, #FALSE,
                                   seed = 35) {
  set.seed(seed)
  threshold_method <- match.arg(threshold_method)
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  site_names <- names(datasources)
  n_sites <- length(site_names)
  
  # Determine if adaptive or fixed for each parameter
  adaptive_theta <- is.null(theta_every)
  adaptive_us <- is.null(us_every)
  adaptive_swap <- is.null(sync_every)
  
  # Set initial values for adaptive, or fixed values for fixed
  current_theta_every <- if (adaptive_theta) initial_every else theta_every
  current_us_every <- if (adaptive_us) initial_every else us_every
  current_sync_every <- if (adaptive_swap) (initial_every * 2) else sync_every
  
  # ---- 1. WARM-START ----
  message("Warm-start initialization (", warmup_burn, " iterations)")
  message("  θ: ", if (adaptive_theta) "ADAPTIVE" else "FIXED (every ", theta_every, ")")
  message("  u: ", if (adaptive_us) "ADAPTIVE" else "FIXED (every ", us_every, ")")
  message("  swap: ", if (adaptive_swap) "ADAPTIVE" else "FIXED (every ", sync_every, ")")
  
  site_meta <- DSI::datashield.aggregate(datasources,
      call("semiOPBARTLocalInitDS", data.name, num_tree, k, 5, state_name, allow_local_ecdf, seed = seed))
  
  p <- site_meta[[1]]$p
  J <- site_meta[[1]]$J
  if (!all(vapply(site_meta, function(m) m$p == p && m$J == J, logical(1))))
    stop("sites disagree on p (ncol W) or J (number of categories)")
  
  theta <- rep(0, p)
  us <- 0:(J - 2); us[1] <- 0
  
  # Run warmup with communication every iteration
  if (warmup_burn > 0) {
    warmup_result <- ds.semiOPBARTWarmupBlock(
      conns = datasources,
      site_names = site_names,
      data.name = data.name,
      num_tree = num_tree,
      k = k,
      num_iterations = warmup_burn,
      rotate_every = 1,  # Communicate every iteration during warmup
      n_swap = n_swap,
      theta = theta,
      us = us,
      p = p,
      J = J,
      nfilter_threshold = nfilter_threshold,
      threshold_method = threshold_method,
      state_name = state_name,
      n_sites = n_sites,
      adapt_rate = 0,
      seed = seed
    )
    
    theta <- warmup_result$theta
    us <- warmup_result$us
    
    message("Warm-up complete. Final theta range: [", 
            round(min(theta), 3), ", ", round(max(theta), 3), "]")
  }
  
  # ---- 2. MAIN LOOP ----
  message("Main MCMC")
  
  N <- num_burn + num_save
  theta_draws <- matrix(NA_real_, N, p)
  us_draws <- matrix(NA_real_, N, length(us))
  
  total_sweeps_done <- 0
  
  # Track theta changes for adaptive scheduling
  all_theta_changes <- numeric()
  all_us_changes <- numeric()
  
  # Track when each parameter was last updated
  last_theta_sync <- 0
  last_us_sync <- 0
  last_swap <- 0
  
while (total_sweeps_done < N) {
    
    # Determine if this block should update each parameter
    sweeps_remaining <- N - total_sweeps_done
    
    # Calculate how many sweeps until next required update
    next_theta_sync <- last_theta_sync + current_theta_every
    next_us_sync <- last_us_sync + current_us_every
    next_swap <- last_swap + current_sync_every
    
    # Number of sweeps to run before the next scheduled communication
    next_event <- min(next_theta_sync, next_us_sync, next_swap, N)
    num_local_sweeps <- min(next_event - total_sweeps_done, sweeps_remaining)
    
    # ---- FIX: Ensure num_local_sweeps is an integer ----
    num_local_sweeps <- max(1, as.integer(round(num_local_sweeps)))
    
    if (num_local_sweeps <= 0) {
      num_local_sweeps <- 1
    }
    
    # Determine which parameters to update in this block
    update_theta <- (total_sweeps_done + num_local_sweeps >= next_theta_sync)
    update_us <- (total_sweeps_done + num_local_sweeps >= next_us_sync)
    perform_swap <- (total_sweeps_done + num_local_sweeps >= next_swap && n_sites > 1)
    
    # ---- 2a. Run local MCMC on all sites ----
    site_results <- DSI::datashield.aggregate(datasources,
        call("semiOPBARTLocalMCMCImproved",
             num_local_sweeps,
             semiOPBART_toSerialize(theta),
             semiOPBART_toSerialize(us),
             update_theta,
             update_us,
             perform_swap,
             n_swap,
             "null",
             "null",
             nfilter_threshold,
             J,
             threshold_method,
             state_name,
             adapt_rate,
             seed = seed))
    
    total_sweeps_done <- total_sweeps_done + num_local_sweeps
    
    # ---- Extract theta changes for adaptive scheduling ----
    theta_changes <- sapply(site_results, function(res) {
      if (is.list(res) && !is.null(res$theta_change) && is.finite(res$theta_change)) {
        return(res$theta_change)
      } else {
        return(NA_real_)
      }
    })
    
    valid_changes <- theta_changes[!is.na(theta_changes) & is.finite(theta_changes)]
    if (length(valid_changes) > 0) {
      all_theta_changes <- c(all_theta_changes, valid_changes)
      if (length(all_theta_changes) > 100) {
        all_theta_changes <- tail(all_theta_changes, 100)
      }
    }
    
    # ---- 2b. AGGREGATE THETA (if scheduled) ----
    if (update_theta) {
      WtW_list <- lapply(site_results, `[[`, "WtW")
      WtZr_list <- lapply(site_results, `[[`, "WtZr")
      valid_idx <- which(!sapply(WtW_list, is.null))
      
      if (length(valid_idx) > 0) {
        WtW <- Reduce(`+`, WtW_list[valid_idx])
        WtZr <- Reduce(`+`, WtZr_list[valid_idx])
        theta_hat <- solve(WtW) %*% WtZr
        theta_sigma <- solve(WtW)
        theta <- as.numeric(mvtnorm::rmvnorm(1, theta_hat, theta_sigma))
        last_theta_sync <- total_sweeps_done
      }
      
      # Update theta communication schedule
      if (adaptive_theta && length(all_theta_changes) > 0) {
        current_theta_every <- adaptive_scheduler(
          theta_changes = all_theta_changes,
          current_every = current_theta_every,
          min_every = min_every,
          max_every = max_every
        )
      }
    }
    
    # ---- 2c. AGGREGATE THRESHOLDS (if scheduled) ----
    if (update_us) {
      th_stats <- lapply(site_results, `[[`, "th_stats")
      th_stats <- th_stats[!sapply(th_stats, is.null)]
      if (length(th_stats) > 0) {
        us <- update_thresholds_from_site_stats(th_stats, us)
        last_us_sync <- total_sweeps_done
      }
      
      # Update us communication schedule
      if (adaptive_us && length(all_theta_changes) > 0) {
        current_us_every <- adaptive_scheduler(
          theta_changes = all_theta_changes,
          current_every = current_us_every,
          min_every = min_every,
          max_every = max_every
        )
      }
    }
    
    # ---- 2d. TREE SWAPPING (if scheduled) ----
    if (perform_swap && n_sites > 1) {
      has_trees <- sapply(site_results, function(res) {
        !is.null(res$trees_Serialize) && res$trees_Serialize != "null"
      })
      
      if (sum(has_trees) > 1) {
        valid_sites <- names(site_results)[has_trees]
        export_idx <- setNames(
          lapply(valid_sites, function(s) sort(sample(0:(num_tree - 1), n_swap))),
          valid_sites)
        
        export_trees <- setNames(lapply(valid_sites, function(s) {
          all_trees <- semiOPBART_fromSerialize(site_results[[s]]$trees_Serialize)
          all_trees[export_idx[[s]] + 1]
        }), valid_sites)
        
        # FIX: was `for (k in seq_along(valid_sites))`, silently shadowing
        # this function's own `k` parameter (BART shrinkage) -- see the
        # identical fix in ds.semiOPBARTWarmupBlock() and ds.semiOPBART().
        for (swap_i in seq_along(valid_sites)) {
          sender <- valid_sites[swap_i]
          receiver <- valid_sites[if (swap_i == length(valid_sites)) 1 else swap_i + 1]
          
          DSI::datashield.aggregate(datasources[receiver],
              call("semiOPBARTLocalSwapDS", 
                   semiOPBART_toSerialize(export_idx[[receiver]]),
                   semiOPBART_toSerialize(export_trees[[sender]]),
                   state_name,
                   seed = seed  ))
        }
        last_swap <- total_sweeps_done
      }
      
      # Update swap schedule
      if (adaptive_swap) {
        # Swap less frequently over time when adaptive
        current_sync_every <- min(max_every, current_sync_every * 1.05)
      }
    }
    
    # ---- 2e. Store draws ----
    block_start <- total_sweeps_done - num_local_sweeps + 1
    block_end <- total_sweeps_done
    
    for (idx in block_start:block_end) {
      theta_draws[idx, ] <- theta
      us_draws[idx, ] <- us
    }
    
    # Log progress
    if (total_sweeps_done %% 100 == 0 || total_sweeps_done == N) {
      msg <- paste0("Progress: ", total_sweeps_done, "/", N)
      if (adaptive_theta) {
        msg <- paste0(msg, ", θ_every=", current_theta_every)
      }
      if (adaptive_us) {
        msg <- paste0(msg, ", u_every=", current_us_every)
      }
      if (adaptive_swap) {
        msg <- paste0(msg, ", sync_every=", current_sync_every)
      }
      message(msg)
    }
  }
  
  keep <- (num_burn + 1):N
  list(theta = theta_draws[keep, , drop = FALSE],
       us = us_draws[keep, , drop = FALSE],
       data.name = data.name,
       num_tree = num_tree, 
       k = k,
       site_names = site_names, 
       state_name = state_name,
       adaptive_theta = adaptive_theta,
       adaptive_us = adaptive_us,
       adaptive_swap = adaptive_swap,
       final_theta_every = if (adaptive_theta) current_theta_every else theta_every,
       final_us_every = if (adaptive_us) current_us_every else us_every,
       final_sync_every = if (adaptive_swap) current_sync_every else sync_every)
}




#' ds.semiOPBARTTrainEImproved - Architecture E with optional adaptive communication
#' 
#' Architecture E: Only theta/us are synchronized periodically.
#' Trees remain entirely local at each site.
#'
#' @param data.name  TRAIN object name
#' @param datasources list of DataSHIELD connections
#' @param num_tree number of trees per site
#' @param k BART shrinkage parameter
#' @param num_burn number of burn-in iterations
#' @param num_save number of saving iterations
#' @param sync_every frequency for theta/us aggregation. If NULL (default), 
#'                   uses adaptive communication. If numeric, uses fixed 
#'                   periodic communication.
#' @param initial_every initial communication frequency (adaptive only)
#' @param min_every minimum communication interval (adaptive only)
#' @param max_every maximum communication interval (adaptive only)
#' @param warmup_burn number of warmup iterations
#' @param adapt_rate local adaptation rate during long gaps
#' @param threshold_method "percentile" or "exact"
#' @param state_name where to find/save local state
#' @export
ds.semiOPBARTTrainEImproved <- function(data.name = "D_train",
                                         datasources = NULL,
                                         num_tree = 20,
                                         k = 1,
                                         num_burn = 1000,
                                         num_save = 1000,
                                         sync_every = NULL,
                                         initial_every = 50,
                                         min_every = 5,
                                         max_every = 200,
                                         warmup_burn = 100,
                                         adapt_rate = 0.1,
                                         nfilter_threshold = 5,
                                         threshold_method = c("percentile", "exact"),
                                         state_name = ".semiOPBART_state",
                                         seed = 35) {
  set.seed(seed)
  threshold_method <- match.arg(threshold_method)
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  site_names <- names(datasources)
  n_sites <- length(site_names)
  
  # Determine if adaptive or fixed
  adaptive_sync <- is.null(sync_every)
  current_sync_every <- if (adaptive_sync) initial_every else sync_every
  
  # ---- 1. WARM-START ----
  message("Warm-start initialization (", warmup_burn, " iterations)")
  message("  sync: ", if (adaptive_sync) "ADAPTIVE" else "FIXED (every ", sync_every, ")")
  
  site_meta <- DSI::datashield.aggregate(datasources,
      call("semiOPBARTLocalInitEDS", data.name, num_tree, k, 5, state_name, seed = seed))
  p <- site_meta[[1]]$p
  J <- site_meta[[1]]$J
  if (!all(vapply(site_meta, function(m) m$p == p && m$J == J, logical(1))))
    stop("sites disagree on p (ncol W) or J (number of categories)")
  
  theta <- rep(0, p)
  us <- 0:(J - 2); us[1] <- 0
  
  if (warmup_burn > 0) {
    warmup_result <- ds.semiOPBARTWarmupBlockE(
      conns = datasources,
      site_names = site_names,
      data.name = data.name,
      num_tree = num_tree,
      k = k,
      num_iterations = warmup_burn,
      sync_every = 1,
      theta = theta,
      us = us,
      p = p,
      J = J,
      nfilter_threshold = nfilter_threshold,
      threshold_method = threshold_method,
      state_name = state_name,
      adapt_rate = 0,
      seed = seed
    )
    theta <- warmup_result$theta
    us <- warmup_result$us
  }
  
  # ---- 2. MAIN LOOP ----
  message("Main MCMC")
  
  N <- num_burn + num_save
  theta_draws <- matrix(NA_real_, N, p)
  us_draws <- matrix(NA_real_, N, length(us))
  
  total_sweeps_done <- 0
  all_theta_changes <- numeric()
  
  while (total_sweeps_done < N) {
    
    sweeps_remaining <- N - total_sweeps_done
    num_local_sweeps <- min(current_sync_every, sweeps_remaining)
    
    # ---- 2a. Run local MCMC ----
    # Architecture E: No tree swapping ever
    site_results <- DSI::datashield.aggregate(datasources,
        call("semiOPBARTLocalMCMCImproved",
             num_local_sweeps,
             semiOPBART_toSerialize(theta),
             semiOPBART_toSerialize(us),
             TRUE,  # update_theta
             TRUE,  # update_us
             FALSE, # no tree swapping
             0,
             "null",
             "null",
             nfilter_threshold,
             J,
             threshold_method,
             state_name,
             adapt_rate,
             seed = seed))
    
    total_sweeps_done <- total_sweeps_done + num_local_sweeps
    
    # ---- Extract theta changes for adaptive scheduling ----
    theta_changes <- sapply(site_results, function(res) {
      if (is.list(res) && !is.null(res$theta_change) && is.finite(res$theta_change)) {
        return(res$theta_change)
      } else {
        return(NA_real_)
      }
    })
    
    valid_changes <- theta_changes[!is.na(theta_changes) & is.finite(theta_changes)]
    if (length(valid_changes) > 0) {
      all_theta_changes <- c(all_theta_changes, valid_changes)
      if (length(all_theta_changes) > 100) {
        all_theta_changes <- tail(all_theta_changes, 100)
      }
    }
    
    # ---- 2b. AGGREGATE THETA ----
    WtW_list <- lapply(site_results, `[[`, "WtW")
    WtZr_list <- lapply(site_results, `[[`, "WtZr")
    valid_idx <- which(!sapply(WtW_list, is.null))
    
    if (length(valid_idx) > 0) {
      WtW <- Reduce(`+`, WtW_list[valid_idx])
      WtZr <- Reduce(`+`, WtZr_list[valid_idx])
      theta_hat <- solve(WtW) %*% WtZr
      theta_sigma <- solve(WtW)
      theta <- as.numeric(mvtnorm::rmvnorm(1, theta_hat, theta_sigma))
    }
    
    # ---- 2c. AGGREGATE THRESHOLDS ----
    th_stats <- lapply(site_results, `[[`, "th_stats")
    th_stats <- th_stats[!sapply(th_stats, is.null)]
    if (length(th_stats) > 0) {
      us <- update_thresholds_from_site_stats(th_stats, us)
    }
    
    # ---- 2d. Update communication schedule ----
    if (adaptive_sync && length(all_theta_changes) > 0) {
      current_sync_every <- adaptive_scheduler(
        theta_changes = all_theta_changes,
        current_every = current_sync_every,
        min_every = min_every,
        max_every = max_every
      )
    }
    
    # ---- 2e. Store draws ----
    block_start <- total_sweeps_done - num_local_sweeps + 1
    block_end <- total_sweeps_done
    
    for (idx in block_start:block_end) {
      theta_draws[idx, ] <- theta
      us_draws[idx, ] <- us
    }
    
    if (total_sweeps_done %% 100 == 0 || total_sweeps_done == N) {
      msg <- paste0("Progress: ", total_sweeps_done, "/", N)
      if (adaptive_sync) {
        msg <- paste0(msg, ", sync_every=", current_sync_every)
      }
      message(msg)
    }
  }
  
  keep <- (num_burn + 1):N
  list(theta_draws = theta_draws[keep, , drop = FALSE],
       us_draws = us_draws[keep, , drop = FALSE],
       theta_mean = colMeans(theta_draws[keep, , drop = FALSE]),
       us_mean = colMeans(us_draws[keep, , drop = FALSE]),
       site_names = site_names,
       data.name = data.name,
       num_tree = num_tree,
       k = k,
       state_name = state_name,
       adaptive_sync = adaptive_sync,
       final_sync_every = if (adaptive_sync) current_sync_every else sync_every)
}


#' Warmup block for Architecture D - FIXED
ds.semiOPBARTWarmupBlock <- function(conns, site_names, data.name, num_tree, k,
                                      num_iterations, rotate_every, n_swap,
                                      theta, us, p, J, nfilter_threshold,
                                      threshold_method, state_name, n_sites,
                                      adapt_rate = 0, seed = 35) {
  
  total_done <- 0
  set.seed(seed)
  while (total_done < num_iterations) {
    sweeps_remaining <- num_iterations - total_done
    num_local_sweeps <- min(rotate_every, sweeps_remaining)
    
    # ---- FIX: Ensure num_local_sweeps is an integer ----
    num_local_sweeps <- as.integer(round(num_local_sweeps))
    
    is_last_block <- (total_done + num_local_sweeps >= num_iterations)
    
    perform_swap <- (n_sites > 1 && 
                     num_local_sweeps == rotate_every && 
                     !is_last_block)
    
    site_results <- DSI::datashield.aggregate(conns,
        call("semiOPBARTLocalMCMCImproved",
             num_local_sweeps,
             semiOPBART_toSerialize(theta),
             semiOPBART_toSerialize(us),
             TRUE,  # update_theta
             TRUE,  # update_us
             perform_swap,
             n_swap,
             "null",
             "null",
             nfilter_threshold,
             J,
             threshold_method,
             state_name,
             adapt_rate,
             seed = seed))
    
    total_done <- total_done + num_local_sweeps
    
    # Aggregate theta - with NULL checks
    WtW_list <- lapply(site_results, `[[`, "WtW")
    WtZr_list <- lapply(site_results, `[[`, "WtZr")
    valid_idx <- which(!sapply(WtW_list, is.null))
    
    if (length(valid_idx) > 0) {
      WtW <- Reduce(`+`, WtW_list[valid_idx])
      WtZr <- Reduce(`+`, WtZr_list[valid_idx])
      theta_hat <- solve(WtW) %*% WtZr
      theta_sigma <- solve(WtW)
      theta <- as.numeric(mvtnorm::rmvnorm(1, theta_hat, theta_sigma))
    }
    
    # Aggregate thresholds
    th_stats <- lapply(site_results, `[[`, "th_stats")
    th_stats <- th_stats[!sapply(th_stats, is.null)]
    if (length(th_stats) > 0) {
      us <- update_thresholds_from_site_stats(th_stats, us)
    }
    
    # Tree swapping
    if (perform_swap && n_sites > 1) {
      has_trees <- sapply(site_results, function(res) {
        !is.null(res$trees_Serialize) && res$trees_Serialize != "null"
      })
      
      if (sum(has_trees) > 1) {
        valid_sites <- names(site_results)[has_trees]
        export_idx <- setNames(
          lapply(valid_sites, function(s) sort(sample(0:(num_tree - 1), n_swap))),
          valid_sites)
        
        export_trees <- setNames(lapply(valid_sites, function(s) {
          all_trees <- semiOPBART_fromSerialize(site_results[[s]]$trees_Serialize)
          all_trees[export_idx[[s]] + 1]
        }), valid_sites)
        
        # FIX: was `for (k in seq_along(valid_sites))`, which silently
        # shadowed this function's own `k` parameter (BART shrinkage) for
        # the rest of the function -- harmless only because nothing below
        # re-reads k, but a real latent bug if that ever changes.
        for (swap_i in seq_along(valid_sites)) {
          sender <- valid_sites[swap_i]
          receiver <- valid_sites[if (swap_i == length(valid_sites)) 1 else swap_i + 1]
          DSI::datashield.aggregate(conns[receiver],
              call("semiOPBARTLocalSwapDS", 
                   semiOPBART_toSerialize(export_idx[[receiver]]),
                   semiOPBART_toSerialize(export_trees[[sender]]),
                   state_name,
                   seed = seed))
        }
      }
    }
  }
  
  list(theta = theta, us = us)
}


# dsSemiOPBARTClient.R - Fixed version

#' Adaptive communication scheduler - Fixed
#' 
#' @param theta_changes  numeric vector of theta changes from each site
#' @param current_every  current communication frequency
#' @param min_every      minimum frequency (most communication)
#' @param max_every      maximum frequency (least communication)
#' @param target_change  target relative change threshold
#' @param momentum       smoothing parameter (0-1)
adaptive_scheduler <- function(theta_changes, 
                               current_every,
                               min_every = 5,
                               max_every = 200,
                               target_change = 0.01,
                               momentum = 0.3) {
  
  # Filter out NULL, NA, and non-finite values
  valid_changes <- theta_changes[!is.na(theta_changes) & 
                                  is.finite(theta_changes) & 
                                  !is.null(theta_changes)]
  
  if (length(valid_changes) == 0) {
    # No valid changes - maintain current frequency
    return(current_every)
  }
  
  # Average change across sites
  avg_change <- mean(valid_changes, na.rm = TRUE)
  
  # If change is zero or very small, gradually increase every
  if (avg_change < 1e-10) {
    new_every <- min(max_every, current_every * 1.1)
    return(round(new_every))
  }
  
  # Adjust frequency based on change
  # Higher change -> lower every (more frequent communication)
  ratio <- target_change / avg_change
  new_every <- current_every * ratio
  
  # Apply bounds
  new_every <- max(min_every, min(max_every, new_every))
  
  # Apply momentum for smooth transitions
  new_every <- momentum * new_every + (1 - momentum) * current_every
  
  return(round(new_every))
}

#' Warmup block for Architecture E - FIXED
ds.semiOPBARTWarmupBlockE <- function(conns, site_names, data.name, num_tree, k,
                                       num_iterations, sync_every,
                                       theta, us, p, J, nfilter_threshold,
                                       threshold_method, state_name,
                                       adapt_rate = 0, seed = 35) {
  
    
  total_done <- 0
  set.seed(seed)
  while (total_done < num_iterations) {
    sweeps_remaining <- num_iterations - total_done
    num_local_sweeps <- min(sync_every, sweeps_remaining)
    
    # ---- FIX: Ensure num_local_sweeps is an integer ----
    num_local_sweeps <- as.integer(round(num_local_sweeps))
    
    site_results <- DSI::datashield.aggregate(conns,
        call("semiOPBARTLocalMCMCImproved",
             num_local_sweeps,
             semiOPBART_toSerialize(theta),
             semiOPBART_toSerialize(us),
             TRUE,  # update_theta
             TRUE,  # update_us
             FALSE, # no swapping
             0,
             "null",
             "null",
             nfilter_threshold,
             J,
             threshold_method,
             state_name,
             adapt_rate,
             seed = seed))
    
    total_done <- total_done + num_local_sweeps
    
    # Aggregate theta - with NULL checks
    WtW_list <- lapply(site_results, `[[`, "WtW")
    WtZr_list <- lapply(site_results, `[[`, "WtZr")
    valid_idx <- which(!sapply(WtW_list, is.null))
    
    if (length(valid_idx) > 0) {
      WtW <- Reduce(`+`, WtW_list[valid_idx])
      WtZr <- Reduce(`+`, WtZr_list[valid_idx])
      theta_hat <- solve(WtW) %*% WtZr
      theta_sigma <- solve(WtW)
      theta <- as.numeric(mvtnorm::rmvnorm(1, theta_hat, theta_sigma))
    }
    
    # Aggregate thresholds
    th_stats <- lapply(site_results, `[[`, "th_stats")
    th_stats <- th_stats[!sapply(th_stats, is.null)]
    if (length(th_stats) > 0) {
      us <- update_thresholds_from_site_stats(th_stats, us)
    }
  }
  
  list(theta = theta, us = us)
}







#' @export
ds.semiOPBARTTrain <- ds.semiOPBARTImproved    # explicit train-only alias, see
                                       # dsSemiOPBARTTrainTest.R   

#' @export
ds.semiOPBARTTrainE <- ds.semiOPBARTTrainEImproved   # explicit train-only alias, see
                                       # dsSemiOPBARTTrainTest.R   
