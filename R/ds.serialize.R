# dsSemiOPBARTSerialize.R
# ---------------------------------------------------------------------------
# DataSHIELD's transport only allows STRING arguments to reach the server:
# every datashield.aggregate()/datashield.assign() call is a function name
# plus a set of arguments that must survive being sent to Opal as text and
# reconstructed there. A bare scalar number or short string is fine as a
# literal; anything else -- numeric vectors of length > 1, matrices, lists,
# data frames, nested tree structures -- must be explicitly Serialize-encoded by
# the CLIENT before the call, and Serialize-decoded as the first line of the
# corresponding SERVER function. This file is the single place that logic
# lives; every other file in this package sources it and uses these two
# functions at every such boundary. (Return values coming BACK from a
# server function are not subject to this restriction -- DSI's aggregate
# mechanism serializes results natively -- so decoding is only ever needed
# on the server side of an incoming call, and encoding only ever needed on
# the client side of an outgoing call.)
# ---------------------------------------------------------------------------


#' Ported verbatim from the real smopbart() source (confirmed, not
#' reconstructed): groups a dummyVars object's expanded dummy columns by
#' their ORIGINAL categorical variable, so CSBart's Dirichlet variable-
#' selection prior treats a multi-level categorical as ONE logical
#' variable rather than letting its dummy columns compete independently.
#' Shared here (rather than duplicated) because it's needed anywhere a
#' fresh Hypers object is built from a dummyVars object: Init
#' (dsSemiOPBARTBase.R) and the never-trained-site path
#' (dsSemiOPBARTExternalPredict.R) both need it -- both previously never
#' set hypers$group at all, silently diverging from the real algorithm
#' whenever any x_feature had more than 2 levels.
dummy_assign <- function(dummy) {
  terms <- attr(dummy$terms, "term.labels")
  group <- list()
  j <- 0
  for (k in terms) {
    if (k %in% dummy$facVars) {
      group[[k]] <- rep(j, length(dummy$lvls[[k]]))
    } else {
      group[[k]] <- j
    }
    j <- j + 1
  }
  do.call(c, group)
}

#' Encode any R object (vector, matrix, list, nested list, data.frame) to a
#' Serialize string safe to embed as a DataSHIELD call argument.
#'
#' Matrices get an explicit dim/dimnames wrapper because plain
#' toSerialize()/fromSerialize() on a matrix silently degrades to a list of rows and
#' loses shape -- semiOPBART_fromSerialize() below knows how to reverse this.
#' @export
semiOPBART_toSerialize  <- function(x) {
    paste(format(serialize(x, NULL)), collapse = "")
}

#' Encode any R object as Hex Serialized string safe to embed as a DataSHIELD call argument.

#' Inverse of semiOPBART_toSerialize(); reconstructs matrices from the wrapper.
#' Always call this as the FIRST line of any server function argument that
#' was Serialize-encoded on the client side -- never operate on the raw string.
#' Decode base64(Serialize)
#' @export
semiOPBART_fromSerialize<- function(hex) {
    raw <- as.raw(strtoi(
        substring(hex,
                  seq(1, nchar(hex), 2),
                  seq(2, nchar(hex), 2)),
        16L
    ))
    unserialize(raw)
}



#' Serialize an ecdf() (or the min-max fallback closure used for 2-level
#' columns in dsSemiOPBARTDataPrep.R's local_ecdf path) to Serialize. This is
#' EXACT, not approximate: R's ecdf() is fundamentally a step function --
#' a table of breakpoints and cumulative probabilities -- so shipping
#' knots(f) and f(knots(f)) losslessly reconstructs it. Used only by
#' Option C's opt-in external-prediction path (dsSemiOPBARTLocalCombine.R),
#' where a reference site's own normalization has to travel with its
#' forest for a never-trained site to use either meaningfully.
semiOPBART_ecdfToSerialize <- function(f) {
  x <- tryCatch(stats::knots(f), error = function(e) NULL)
  if (is.null(x)) {
    # not a stepfun-based ecdf -- must be the 2-level min-max closure from
    # dsSemiOPBARTDataPrep.R; probe it at 0/1 to recover a/b algebraically
    # (y = (x-a)/(b-a) => f(0) = -a/(b-a), f(1) = (1-a)/(b-a))
    f0 <- f(0); f1 <- f(1)
    b_minus_a <- 1 / (f1 - f0)
    a <- -f0 * b_minus_a
    return(semiOPBART_toSerialize(list(.minmax = TRUE, a = a, b = a + b_minus_a)))
  }
  semiOPBART_toSerialize(list(.minmax = FALSE, x = x, y = f(x)))
}

#' Reconstruct a normalization closure from semiOPBART_ecdfToSerialize()'s output.
semiOPBART_ecdfFromSerialize <- function(Serialize_str) {
  spec <- semiOPBART_fromSerialize(Serialize_str)
  if (isTRUE(spec$.minmax)) {
    a <- spec$a; b <- spec$b
    return(function(y) (y - a) / (b - a))
  }
  x <- spec$x; y <- spec$y
  function(newx) stats::approx(x, y, xout = newx, method = "constant",
                                rule = 2, ties = "ordered")$y
}

#' Encode a Forest's exported tree list (from the get_trees() Rcpp method
#' added in csbart_forest_serialize_patch.cpp) to Serialize. Kept as a named
#' wrapper rather than a bare call to semiOPBART_toSerialize() so call sites read
#' clearly, and so the encoding strategy for trees specifically can change
#' later (e.g. a more compact array-based format) without touching every
#' caller.
# semiOPBART_treesToSerialize <- function(forest, tree_idx) {
#   semiOPBART_toSerialize(forest$get_trees(tree_idx))
# }

semiOPBART_treesToSerialize <- function(forest) {
  semiOPBART_toSerialize(forest)
}

#' Decode a Serialize tree bundle and load it into `forest` at `tree_idx`
#' (overwriting whatever was previously at those indices).
semiOPBART_treesFromSerialize <- function(forest, Serialize_str, tree_idx) {
  # trees <- Serializelite::fromSerialize(Serialize_str, simplifyVector = FALSE)
  # forest$set_trees(trees, tree_idx)
  trees <- semiOPBART_fromSerialize(Serialize_str)
  forest$set_trees(trees, tree_idx)
  invisible(NULL)
}

semiOPBART_treesFromSerialize_check <- function(Serialize_str) {

  if (
    length(Serialize_str) != 1L ||
    is.na(Serialize_str) ||
    !nzchar(Serialize_str)
  ) {
    stop(
      "semiOPBART_treesFromSerialize(): ",
      "Serialize_str must be one non-empty, non-NA string"
    )
  }

  trees <- semiOPBART_fromSerialize(
    Serialize_str
  )

  if (is.null(trees)) {
    stop(
      "semiOPBART_treesFromSerialize(): ",
      "deserialization returned NULL"
    )
  }

  trees
}