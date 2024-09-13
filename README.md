DNAcycP R package 
================

**Maintainer**: Ji-Ping Wang, \<<jzwang@northwestern.edu>\>; Brody Kendall \<<curtiskendall2025@u.northwestern.edu>\>; Keren Li, \<<keren.li@northwestern.edu>\>

TODO: license

**License**:

**Cite DNAcycP package**:

TODO: Update citation when applicable

Li, K., Carroll, M., Vafabakhsh, R., Wang, X.A. and Wang, J.-P., DNAcycP: A Deep Learning Tool for DNA Cyclizability Prediction, *Nucleic Acids Research*, 2021

## What is DNAcycP?

**DNAcycP**, short for **DNA** **cyc**lizablity **P**rediction, is an R package for accurate prediction of DNA intrinsic cyclizablity score. It was built upon a deep learning architecture with a hybrid of Inception and Residual network structure and an LSTM layer. The original DNAcycP was trained based on loop-seq data from Basu et al 2021 (see below). An updated version (DNAcycP2) was trained based on smoothed predictions of this loop-seq data. The predicted score (for either DNAcycP or DNAcycP2), termed **C-score** achieves high accuracy compared to the experimentally measured cyclizablity score by loop-seq assay.

## Available format of DNAcycP

TODO: update reference to R package
TODO: update web server

DNAcycP is available in three formats: A web server available at http://DNAcycP.stats.northwestern.edu for real-time prediction and visualization of C-score up to 20K bp, a standalone Python package available for free download from https://github.com/jipingw/DNAcycP, and an R package.

## Architecture of DNAcycP

The core of DNAcycP is a deep learning architecture mixed with an Inception-ResNet structure and an LSTM layer (IR+LSTM, Fig 1b) that processes the sequence and its reverse complement separately, the results from which are averaged and detrended to reach the predicted intrinsic score. (Fig 1a).

IR+LSTM starts with a convolutional layer for dimension reduction such that the encoded sequence space is reduced from 2D to 1D. The output is fed into an inception module that contains two parallel branches, each having two sequentially connected convolutional layers with branch-specific kernels to capture sequence features of different scale. The first branch has kernel dimension 3x1 for both layers and the second has kernel dimension 11x1 and 21x1 sequentially. The output of the inception module is combined by concatenation and added back to the input of the inception module to form a short circuit or residual network. Finally, the IR+LSTM concludes with a dense layer to predict output with linear activation function. 

![A diagram of DNAcycP.](./figures/Figure1.png)

## DNAcycP required R packages:

* `basilisk`
* `reticulate`

## Installation

Current best practice is to install via `devtools` and github:

```r
devtools::install_github("brodykendall/dnacycp-R")
```


## Usage

Upon successful installation, the DNAcycP r package supports the input sequence in two formats: FASTA format (with sequence name line beginning with “>”) or directly as an R object. Unlike in the web server version where only one sequence is allowed in input for prediction, the R package allows multiple sequences in the same input file/object. In particular for the R object format, each sequence (which can be of length >= 50bp) in the file is regarded as one input sequence for prediction, however the computation is most efficient when every has length exactly 50bp.

The two main functions in the package are `cycle` and `cycle_fasta`, both of which perform cyclizability. The main difference between the two functions is the input type: `cycle` takes an R object as its input, while `cycle_fasta` takes a file path as its input. Both take an additional argument `smooth` which determines which model to use in making predictions:
* `smooth=TRUE`: DNAcycP2 (the model trained on smoothed data, recommended) 
* `smooth=FALSE`: DNAcycP (the model trained on the original data).

We provide two simple example files with the package to show proper usage:

### Example 1:

```r
ex1_file <- system.file("data", "ex1.fasta", package = "dnacycp")
ex1_smooth <- dnacycp::cycle_fasta(ex1_file,smooth=TRUE)
ex1_original <- dnacycp::cycle_fasta(ex1_file,smooth=FALSE)
```

`cycle_fasta` takes the file path as input (`ex1_file`)

`smooth=TRUE` specifies that DNAcycP2 be used to make predictions

`smooth=FALSE` specifies that DNAcycP be used to make predictions

The output (`ex1_smooth` or `ex1_original`) is a list with names starting with "cycle"

For example, `ex1.fasta` contains two sequences with IDs "1" and "2" respectively.
Therefore, both `ex1_smooth1` and `ex1_original` will be lists of length 2 with names `cycle_1` and `cycle_2` for the first and second sequences respectively.

 Each item in the list (e.g. `ex1_smooth$cycle_1`) is a data.frame object with three columns: `position`, `c_score_norm`, and `c_score_unnorm`. The `c_score_norm` is the predicted C-score from the model trained based on the standardized loop-seq score (in the case of DNAcycP) or the standardized smoothed intrinsic cyclizability estimate (in the case of DNAcycP2) of the tiling library of Basu et al 2021 (i.e. 0 mean unit variance). When predictions are made using the original DNAcycP, the `c_score_unnorm` is the predicted C-score recovered to the original scale of loop-seq score in the tiling library data from Basu et el 2021. When predictions are made using the updated DNAcycP2 (`smooth=TRUE`), the `c_score_unnorm` is the predicted C-score recovered to the scale of standardized raw cyclizability scores of the tiling library data. The standardized loop-seq score provides two advantages. As loop-seq may be subject to a library-specific constant, standardized C-score is defined with a unified baseline as yeast genome (i.e. 0 mean in yeast genome). Secondly, the C-score provides statisitcal significance indicator, i.e. a C-score of 1.96 indicates 97.5% in the distribution.

### Example 2:

```r
ex2_file <- system.file("data", "ex2.txt", package = "dnacycp")
ex2 <- read.csv(ex2_file, header = FALSE)
ex2_smooth <- dnacycp::cycle(ex2$V1, smooth=TRUE)
ex2_original <- dnacycp::cycle(ex2$V1, smooth=FALSE)
```

`cycle` takes the sequences themselves as input, so we first read the file (`ex2_file`) and then provide the sequences as input (`ex2$V1`)

`smooth=TRUE` specifies that DNAcycP2 be used to make predictions

`smooth=FALSE` specifies that DNAcycP be used to make predictions

The output (`ex2_smooth` or `ex2_original`) is a list with names "norm" and "unnorm" corresponding to normalized and unnormalized C-score.

If every sequence has length exactly 50bp (recommended), both items in the list will be vectors of doubles corresponding to the predicted value for the sequence at the relevant index.

Otherwise (if there as at least one sequence with length >50bp), both items in the list will be lists of vectors corresponding to the predicted values for each subsequence of length 50bp at the relevant list index. For example, as `ex2` contains 100 sequences each of length 250bp, `ex2_smooth$norm[[1]]` contains the normalized C-scores for every 50bp subsequence of the first sequence in `ex2` in order. That is, `ex2_smooth$norm[[1]][1]` corresponds to positions 1-50 of the first sequence in `ex2`, `ex2_smooth$norm[[1]][2]` corresponds to positions 2-51 of the first sequence in `ex2`, and so forth.

## Other References

* Basu, A., Bobrovnikov, D.G., Qureshi, Z., Kayikcioglu, T., Ngo, T.T.M., Ranjan, A., Eustermann, S., Cieza, B., Morgan, M.T., Hejna, M. et al. (2021) Measuring DNA mechanics on the genome scale. Nature, 589, 462-467.