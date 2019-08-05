## FIXME Check that these return same value on coerced objects.



#' @rdname counts
#' @name normcounts
#' @importFrom SingleCellExperiment normcounts
#' @usage normcounts(object, ...)
#' @seealso
#' - `scater::normalizeCounts()`.
#' @export
NULL



## FIXME Need to center at unity?
## Updated 2019-08-04.
`normcounts,SingleCellExperiment` <-  # nolint
    function(object) {
        validObject(object)
        ## Return `normcounts` assay if slotted.
        if (isSubset("normcounts", assayNames(object))) {
            return(assay(x = object, i = "normcounts"))
        }
        ## Otherwise, calculate on the fly using size factors.
        counts <- counts(object)
        sizeFactors <- sizeFactors(object)
        if (is.null(sizeFactors)) {
            stop(paste(
                "Size factors are not defined.",
                "Calculate using `estimateSizeFactors()`."
            ))
        }
        t(t(counts) / sizeFactors)
    }



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
