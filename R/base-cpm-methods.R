## FIXME Add cpm support
## FIXME Check that these return same value on coerced objects.
## FIXME Seurat Relative Counts
## cpm = NormalizeData(object, normalization.method = "RC", scale.factor = 1e6)
## cpm = size factor adjusted / 1e6
## FIXME Explain: Counts-per-million. This is the read count for each gene in each cell, divided by the library size of each cell in millions.
## FIXME `scater::calculateCPM()`.



#' @rdname counts
#' @name counts
#' @importFrom SingleCellExperiment cpm
#' @usage cpm(object, ...)
#' @seealso
#' - `SingleCellExperiment::cpm()`.
#' - `edgeR::cpm()`.
#' - `scater::calculateCPM()`.
#' @export
NULL



## FIXME `edgeR::cpm()`.



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
