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
  # cl <- basiliskStart(env2)
  print("Started basilisk environment 1")
  on.exit(basiliskStop(cl))

  preds <- basiliskRun(cl, fun=function(seqs) {
    print("Entered basilisk Run with env1")
    path_to_python <- system.file("python", package = "dnacycpv2")
    print("Trying to load irlstm:")
    irlstm <- system.file("python/irlstm", package = "dnacycpv2")
    print(irlstm)
    X = reticulate::import_from_path("dnacycpv2_python", path = path_to_python)
    # X <- reticulate::import("dnacycpv2_python")
    print("Imported python scripts")
    res = X$cycle(sequences)
    res
  }, seqs=sequences)

  preds
}


