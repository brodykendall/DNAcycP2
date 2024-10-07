#' Predict Cyclizability
#'
#' This predicts cyclizability for a set of sequences.
#'
#' @param sequences A list or vector of sequences
#' @param smooth Whether to predict smoothed C0 (DNAcycP2) or original C0 (DNAcycP)
#' @return A list of (1) predictions on a normalized scale, (2) predictions on an unnormalized scale
#' @export
#' @importFrom reticulate import_from_path
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
#' @examples
#' # Example usage of cycle
#' cycle(sequences, smooth=TRUE) # where sequences is a list/vector of sequences
cycle <- function(sequences, smooth) {
  cl <- basiliskStart(env1)
  on.exit(basiliskStop(cl))
  
  preds <- basiliskRun(cl, fun=function(seqs) {
    path_to_python <- system.file("python", package = "dnacycp")
    if (smooth) {
      print("Predicting smooth C0:")
      irlstm <- system.file("python/irlstm_smooth", package = "dnacycp")
    }
    else {
      print("Predicting original C0:")
      irlstm <- system.file("python/irlstm", package = "dnacycp")
    }
    X <- reticulate::import_from_path("dnacycp_python", path = path_to_python)
    res <- X$cycle(sequences, irlstm)
    res
  }, seqs=sequences)
  
  preds
}



#' Predict Cyclizability
#'
#' This predicts cyclizability for all subsequences of length 50bp from a .fasta input file.
#'
#' @param input_file .fasta input file path
#' @param smooth Whether to predict smoothed C0 (DNAcycP2) or original C0 (DNAcycP)
#' @return A list of predictions for each ID in the .fasta file. 
#' 
#' Each list item has the following columns: position, c_score_norm (predictions 
#' on a normalized scale), and c_score_unnorm (predictions on an unnormalized
#' scale).
#' 
#' Each list item is named "cycle_{id}" corresponding to the fasta id
#' @export
#' @importFrom reticulate import_from_path
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
#' @examples
#' # Example usage of cycle_fasta
#' cycle_fasta("path/to/fasta/file.fasta",smooth=TRUE)
cycle_fasta <- function(file_path, smooth) {
  cl <- basiliskStart(env1)
  on.exit(basiliskStop(cl))
  
  preds <- basiliskRun(cl, fun=function(input_file) {
    path_to_python <- system.file("python", package = "dnacycp")
    if (smooth) {
      print("Predicting smooth C0:")
      irlstm <- system.file("python/irlstm_smooth", package = "dnacycp")
    }
    else {
      print("Predicting original C0:")
      irlstm <- system.file("python/irlstm", package = "dnacycp")
    }
    X <- reticulate::import_from_path("dnacycp_python", path = path_to_python)
    res <- X$cycle_fasta(input_file, irlstm)
    res
  }, input_file=file_path)
  
  preds
}





