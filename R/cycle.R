utils::globalVariables("model_data")

#' Predict Cyclizability
#'
#' This predicts cyclizability for a set of sequences.
#'
#' Optionally, saves output files (use argument 'save_path_prefix')
#'
#' @param sequences A list or vector of sequences
#' @param smooth Whether to predict smoothed C0 (DNAcycP2) or original C0
#' (DNAcycP)
#' @param save_path_prefix Base path for output files. If it is an empty string,
#' the output files will not be saved (default="")
#' @return A list of predictions for each input sequence.
#' @export
#' @importFrom reticulate import_from_path
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
#' @examples
#' # Example usage of cycle
#' cycle(c("ACTGCTAGTCACTGCTAGTCACTGCTAGTCACTGCTAGTCACTGCTAGTC"), smooth=TRUE)
#' # where sequences is a list/vector of sequences
cycle <- function(sequences, smooth, save_path_prefix="") {
    cl <- basiliskStart(env1)
    on.exit(basiliskStop(cl))
    preds <- basiliskRun(cl, fun=function(seqs) {
        path_to_python <- system.file("python", package = "DNAcycP2")
        irlstm <- restore_model_dir(smooth)
        X <- reticulate::import_from_path(
            "dnacycp_python", path = path_to_python)
        if (inherits(sequences, "AAStringSet") |
            inherits(sequences, "DNAStringSet")
        ) { sequences <- as.character(sequences) }
        res <- X$cycle(sequences, irlstm)
        res
    }, seqs=sequences)
    # Only save outfiles if save_path_prefix argument is supplied:
    if (save_path_prefix != "") {
        # Create file names
        if (smooth) {
            outfile_norm <- file(paste0(
                save_path_prefix, "_C0S_norm.txt"), "w")
            outfile_unnorm <- file(paste0(
                save_path_prefix, "_C0S_unnorm.txt"), "w")
        } else {
            outfile_norm <- file(paste0(
                save_path_prefix, "_C0_norm.txt"), "w")
            outfile_unnorm <- file(paste0(
                save_path_prefix, "_C0_unnorm.txt"), "w")
        }
        # First, save normalized outputs:
        for (row in preds[[2]]) {
            if (is.numeric(row) && length(row) == 1) { s <- as.character(row) }
            else {
                s <- paste(row, collapse = " ")
            }
            writeLines(s, outfile_norm)
        }
        close(outfile_norm)
        # Next, save unnormalized outputs:
        for (row in preds[[3]]) {
            if (is.numeric(row) && length(row) == 1) { s <- as.character(row)}
            else { s <- paste(row, collapse = " ") }
            writeLines(s, outfile_unnorm)
        }
        close(outfile_unnorm)
    }
    preds[[1]]
}



#' Predict Cyclizability
#'
#' This predicts cyclizability for all subsequences of length 50bp from a
#' .fasta input file.
#'
#' Optionally, saves output files (use argument 'save_path_prefix')
#'
#' @param file_path .fasta input file path
#' @param smooth Whether to predict smoothed C0 (DNAcycP2) or original C0
#' (DNAcycP)
#' @param n_cores Number of cores to use for parallel processing (default=1)
#' @param chunk_length Length of sequence that each core will predict on at a
#' given time.
#' (default=100000)
#' @param save_path_prefix Base path for output files. If it is an empty string,
#' the output files will not be saved (default="")
#' @return A list of predictions for each ID in the .fasta file.
#'
#' Each list item has the following columns: position, c_score_norm (
#' predictions on a normalized scale), and c_score_unnorm (predictions on an
#' unnormalized scale).
#'
#' Each list item is named "cycle_$id$" corresponding to the fasta id
#' @export
#' @importFrom reticulate import_from_path
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
#' @importFrom utils write.csv
#'
#' @examples
#' # Create a temporary file
#' temp_file <- tempfile(fileext = ".fasta")
#' writeLines(">1", temp_file)
#' writeLines("ACTGCTAGTCACTGCTAGTCACTGCTAGTCACTGCTAGTCACTGCTAGTC", temp_file)
#'
#' # Example usage of cycle_fasta
#' cycle_fasta(temp_file, smooth=TRUE)
#'
#' # Cleanup
#' unlink(temp_file)
cycle_fasta <- function(file_path, smooth, n_cores=1, chunk_length=100000,
                        save_path_prefix="") {
    cl <- basiliskStart(env1)
    on.exit(basiliskStop(cl))

    preds <- basiliskRun(cl, fun=function(input_file) {
        path_to_python <- system.file("python", package = "DNAcycP2")
        irlstm <- restore_model_dir(smooth)
        X <- reticulate::import_from_path(
            "dnacycp_python", path = path_to_python
        )
        res <- X$cycle_fasta(
            input_file, irlstm, num_threads=as.integer(n_cores),
            chunk_size=as.integer(chunk_length)
        )
        res
    }, input_file=file_path)

    if (save_path_prefix != "") {
        for (i in seq_along(preds)) {
            outfile<-paste0(c(
                paste(c(save_path_prefix,names(preds)[i]),collapse="_"),
                ".txt"), collapse=""
            )
            write.csv(preds[[i]],outfile,row.names=FALSE)
        }
    }

    preds
}



restore_model_dir <- function(smooth) {
    if (smooth) {
        load(system.file("extdata/irlstm_smooth.rda", package = "DNAcycP2"))
        model_dir <- file.path(tempdir(), "irlstm_smooth")
    }
    else {
        load(system.file("extdata/irlstm.rda", package = "DNAcycP2"))
        model_dir <- file.path(tempdir(), "irlstm")
    }

    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

    # Reconstruct files from raw data
    for (rel_path in names(model_data)) {
        out_file <- file.path(model_dir, rel_path)
        dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
        writeBin(model_data[[rel_path]], out_file)
    }

    return(model_dir)
}
