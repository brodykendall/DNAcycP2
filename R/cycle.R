#' Predict Cyclizability
#'
#' This predicts cyclizability for a set of sequences.
#'
#' @param sequences A list or vector of sequences
#' @return A list of (1) predictions on a normalized scale, (2) predictions on an unnormalized scale
#' @export
#' @importFrom reticulate import_from_path
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
cycle <- function(sequences) {
  cl <- basiliskStart(env1)
  on.exit(basiliskStop(cl))

  preds <- basiliskRun(cl, fun=function(seqs) {
    path_to_python <- system.file("python", package = "dnacycpv2")
    irlstm <- system.file("python/irlstm", package = "dnacycpv2")
    X = reticulate::import_from_path("dnacycpv2_python", path = path_to_python)
    res = X$cycle(sequences, irlstm)
    res
  }, seqs=sequences)

  preds
}


