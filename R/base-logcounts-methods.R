## FIXME Check that these return same value on coerced objects.



#' @rdname counts
#' @name logcounts
#' @importFrom SingleCellExperiment logcounts
#' @usage logcounts(object, ...)
#' @export
NULL



## FIXME Consistently return log2 with pseudocount of 1.

## nolint start
## > Seurat:::as.SingleCellExperiment.Seurat
## > object[["RNA"]]@data
## nolint end

## Updated 2019-08-04.
`logcounts,Seurat` <-  # nolint
    function(object, assay = NULL) {
        Seurat::GetAssayData(
            object = object,
            assay = assay,
            slot = "data"
        )
    }



#' @rdname counts
#' @export
setMethod(
    f = "logcounts",
    signature = signature("Seurat"),
    definition = `logcounts,Seurat`
)



`logcounts,cell_data_set` <-  # nolint
    function(object) {
        monocle3::normalized_counts(
            cds = object,
            norm_method = "log",
            pseudocount = 1L
        )
    }



#' @rdname counts
#' @export
setMethod(
    f = "logcounts",
    signature = signature("cell_data_set"),
    definition = `logcounts,cell_data_set`
)
