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



#' Predict Cyclizability
#'
#' This predicts cyclizability for all subsequences of length 50bp from a .fasta input file.
#'
#' @param input_file .fasta input file path
#' @return A list of predictions for each ID in the .fasta file. Each list item has
#'  the following columns: position, c_score_norm (predictions on a normalized scale),
#' and c_score_unnorm (predictions on an unnormalized scale)
#' @export
#' @importFrom reticulate import_from_path
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
cycle_fasta <- function(file_path) {
  cl <- basiliskStart(env1)
  on.exit(basiliskStop(cl))
  
  preds <- basiliskRun(cl, fun=function(input_file) {
    path_to_python <- system.file("python", package = "dnacycpv2")
    irlstm <- system.file("python/irlstm", package = "dnacycpv2")
    X = reticulate::import_from_path("dnacycpv2_python", path = path_to_python)
    res = X$cycle_fasta(input_file, irlstm)
    res
  }, input_file=file_path)
  
  preds
}





