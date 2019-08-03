#' @name counts
#' @inherit SingleCellExperiment::counts
#' @param ... Additional arguments.
NULL



#' @rdname counts
#' @name counts
#' @importFrom SingleCellExperiment counts
#' @usage counts(object, ...)
#' @export
NULL

#' @rdname counts
#' @name logcounts
#' @importFrom SingleCellExperiment logcounts
#' @usage logcounts(object, ...)
#' @export
NULL

#' @rdname counts
#' @name normcounts
#' @importFrom SingleCellExperiment normcounts
#' @usage normcounts(object, ...)
#' @export
NULL



## nolint start
## > Seurat:::as.SingleCellExperiment.Seurat
## > object[["RNA"]]@counts
## nolint end

## Updated 2019-08-03.
`counts,Seurat` <-  # nolint
    function(object, assay = "RNA") {
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



## nolint start
## > Seurat:::as.SingleCellExperiment.Seurat
## > object[["RNA"]]@data
## nolint end

## Updated 2019-08-03.
`logcounts,Seurat` <-  # nolint
    function(object, assay = "RNA") {
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



## FIXME
## This should return size factor adjusted counts, not Seurat's scaled counts
## (i.e. `scale.data`) here.
## Updated 2019-08-03
`normcounts,Seurat` <-  # nolint
    function(object, assay = "RNA") {
        stop("Not yet supported")
    }



#' @rdname counts
#' @export
setMethod(
    f = "normcounts",
    signature = signature("Seurat"),
    definition = `normcounts,Seurat`
)



`normcounts,cell_data_set` <-  # nolint
    function(object) {
        monocle3::normalized_counts(
            cds = object,
            norm_method = "size_only",
            pseudocount = NULL
        )
    }



#' @rdname counts
#' @export
setMethod(
    f = "normcounts",
    signature = signature("cell_data_set"),
    definition = `normcounts,cell_data_set`
)
