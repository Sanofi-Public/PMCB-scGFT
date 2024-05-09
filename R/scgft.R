# Genix package documentation and import directives

#' The scgft package
#'
#' \code{scgft} offers a streamlined computational framework, leveraging
#' the principles of Fourier Transform, to generate synthetic scRNA-seq data in-silico.
#'
#' @section  Synthesis and evaluation of synthetic cells:
#' \itemize{
#'   \item  \link{RunScGFT}:     Generates synthetic scRNA-seq data in-silico
#'   \item  \link{statsScGFT}:   Calculates the likelihood that synthesized
#'   cells will have the same identity as their original counterparts.
#' }
#'
#' @name     scgft
#' @docType  package
#' @references
#' \enumerate{
#'   \item  Nima Nouri: scGFT - single-cell RNA-seq data augmentation
#'   using generative Fourier transformers
#' }
#'
#' @import      methods
#' @importFrom  dplyr           n filter select mutate summarize rename %>%
#' @importFrom  tibble          rownames_to_column
#' @importFrom  Matrix          rowSums colSums Matrix
#' @importFrom  rlang           sym syms
#' @importFrom  stats           fft kmeans sd rnorm
#' @importFrom  stringr         str_which
#' @importFrom  progress        progress_bar
NULL

# Package loading actions
.onAttach <- function(libname, pkgname) {
    msg <- paste("scgft package.",
                 "License and Copyrights in file LICENSE and COPYRIGHTS.",
                 sep="\n")
    packageStartupMessage(msg)
}
