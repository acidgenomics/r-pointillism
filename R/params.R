#' Parameters
#'
#' @name params
#' @keywords internal
#'
#' @param alpha `numeric(1)`.
#'   Adjusted P value ("alpha") cutoff. If left `NULL`, will use the cutoff
#'   defined in the object.
#' @param reducedDim `character(1)`.
#'   Name of reduced dimension matrix slotted in
#'   [`reducedDims()`][SingleCellExperiment::reducedDims].
#'   Includes TNSE, UMAP, PCA, for example.
#'
#' @return No value.
NULL
