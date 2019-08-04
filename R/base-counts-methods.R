## FIXME Add cpm support

## FIXME Seurat Relative Counts
## cpm = NormalizeData(, normalization.method = "RC", scale.factor = 1e6)

## FIXME Add unit tests to check the return on these for coerced data sets.

## cpm = size factor adjusted / 1e6



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

#' @rdname counts
#' @name counts
#' @importFrom SingleCellExperiment cpm
#' @usage cpm(object, ...)
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



## FIXME Confirm this is correct.
## Updated 2019-08-04.
`cpm,Seurat` <-  # nolint
    function(object, assay = NULL) {
        NormalizeData(
            object = object,
            assay = assay,
            ## Relative counts.
            normalization.method = "RC",
            scale.factor = 1e6L
        )
    }



#' @rdname counts
#' @export
setMethod(
    f = "cpm",
    signature = signature("Seurat"),
    definition = `cpm,Seurat`
)



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



## FIXME
## This should return size factor adjusted counts, not Seurat's scaled counts
## (i.e. `scale.data`) here.
## Updated 2019-08-03
`normcounts,Seurat` <-  # nolint
    function(object, assay = NULL) {
        stop("Not yet supported")
        counts <- counts(object, assay = assay)
        ## Apply size factors.
        ## FIXME Rework this step.
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
