#' @name cpm
#' @inherit bioverbs::cpm
#'
#' @seealso
#' - `SingleCellExperiment::cpm()`.
#' - `edgeR::cpm()`.
#' - `scater::calculateCPM()`.
#'
#' @examples
#' data(
#'     SingleCellExperiment,
#'     Seurat,
#'     package = "acidtest"
#' )
#'
#' ## SingleCellExperiment ====
#' object <- SingleCellExperiment
#' object <- estimateSizeFactors(object)
#' x <- cpm(object)
#' class(x)
#' mean(x)
#'
#' ## Seurat ====
#' object <- Seurat
#' cpm(object)
NULL



#' @rdname cpm
#' @name cpm
#' @importFrom bioverbs cpm
#' @usage cpm(object, ...)
#' @export
NULL



## Updated 2019-08-05.
`cpm,SingleCellExperiment` <-  # nolint
    function(object) {
        ## Early return if cpm assay is defined.
        if (isSubset("cpm", assayNames(object))) {
            return(assay(x = object, i = "cpm"))
        }

        ## Otherwise, calculate on the fly.
        assert(is.numeric(sizeFactors(object)))
        calculateCPM(
            object = object,
            exprs_values = "counts",
            use_size_factors = TRUE
        )
    }



#' @rdname cpm
#' @export
setMethod(
    f = "cpm",
    signature = signature("SingleCellExperiment"),
    definition = `cpm,SingleCellExperiment`
)



## Updated 2019-08-05.
`cpm,Seurat` <-  # nolint
    function(object, assay = NULL) {
        object <- NormalizeData(
            object = object,
            assay = assay,
            normalization.method = "RC",
            scale.factor = 1e6L,
            verbose = TRUE
        )
        GetAssayData(
            object = object,
            slot = "data",
            assay = assay
        )
    }



#' @rdname cpm
#' @export
setMethod(
    f = "cpm",
    signature = signature("Seurat"),
    definition = `cpm,Seurat`
)
