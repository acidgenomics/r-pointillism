#' @name counts
#' @inherit SingleCellExperiment::counts
#' @inheritParams Seurat::GetAssayData
#' @param ... Additional arguments.
#'
#' @seealso
#' - `Seurat::GetAssayData()`.
#' - `Seurat::NormalizeData()`.
#' - `Seurat::ScaleData()`.
NULL



#' @rdname counts
#' @name counts
#' @importFrom SingleCellExperiment counts
#' @usage counts(object, ...)
#' @export
NULL



## nolint start
## > Seurat:::as.SingleCellExperiment.Seurat
## > object[["RNA"]]@counts
## nolint end

## Updated 2019-08-04.
`counts,Seurat` <-  # nolint
    function(object, assay = NULL) {
        GetAssayData(
            object = object,
            assay = assay,
            slot = "counts"
        )
    }



#' @rdname counts
#' @export
setMethod(
    f = "counts",
    signature = signature("Seurat"),
    definition = `counts,Seurat`
)
