#' Predict Cyclizability
#'
#' This predicts cyclizability for a set of sequences.
#'
#' @param sequences A list or vector of sequences
#' @return A list of (1) predictions on a normalized scale, (2) predictions on an unnormalized scale
#' @export
#' @importFrom reticulate import
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
cycle <- function(sequences) {
  cl <- basiliskStart(env1)
  on.exit(basiliskStop(cl))

  preds <- basiliskRun(cl, fun=function(seqs) {
    X <- reticulate::import("dnacycpv2")
    res = X$cycle(sequences)
    res
  }, seqs=sequences)

  preds
}


