#' Parameters
#'
#' @name params
#' @keywords internal
#'
#' @param alpha `numeric(1)`.
#'   Adjusted P value ("alpha") cutoff. If left `NULL`, will use the cutoff
#'   defined in the object.
#' @param assay `character(1)`.
#'   Assay name.
#'   For `Seurat`, leave `NULL` or use `"RNA"` to get scRNA-seq data.
#' @param reduction `character(1)`.
#'   Name of dimension reduction matrix slotted in
#'   [`reducedDims()`][SingleCellExperiment::reducedDims].
#'   Includes UMAP, TSNE, PCA, typically, for example.
#'
#' @return No value.
NULL
