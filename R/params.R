#' Parameters
#'
#' @name params
#' @keywords internal
#'
#' @param alpha `numeric(1)`.
#'   Adjusted P value ("alpha") cutoff. If left `NULL`, will use the cutoff
#'   defined in the object.
#' @param reduction `character(1)`.
#'   Name of dimension reduction matrix slotted in
#'   [`reducedDims()`][SingleCellExperiment::reducedDims].
#'   Includes UMAP, TSNE, PCA, typically, for example.
#'
#' @return No value.
NULL
