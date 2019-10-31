#' @name cpm
#' @inherit bioverbs::cpm
#' @keywords internal
#' @note Updated 2019-10-30.
#'
#' @inheritParams acidroxygen::params
#'
#' @seealso
#' - [SingleCellExperiment::cpm()].
#' - [edgeR::cpm()].
#' - [scater::calculateCPM()].
#'
#' @examples
#' data(
#'     Seurat,
#'     SingleCellExperiment,
#'     cell_data_set,
#'     package = "acidtest"
#' )
#'
#' ## SingleCellExperiment ====
#' object <- SingleCellExperiment
#' object <- estimateSizeFactors(object)
#' cpm <- cpm(object)
#' class(cpm)
#' mean(cpm)
#'
#' ## Seurat ====
#' object <- Seurat
#' cpm <- cpm(object)
#' class(cpm)
#' mean(cpm)
#'
#' ## cell_data_set ====
#' ## > object <- cell_data_set
#' ## > cpm <- cpm(object)
#' ## > class(cpm)
#' ## > mean(cpm)
NULL



#' @rdname cpm
#' @name cpm
#' @importFrom bioverbs cpm
#' @usage cpm(object, ...)
#' @export
NULL



## Updated 2019-10-30.
`cpm,SingleCellExperiment` <-  # nolint
    function(object, verbose = FALSE) {
        assert(isFlag(verbose))
        ## Early return if cpm assay is defined.
        if (isSubset("cpm", assayNames(object))) {
            if (isTRUE(verbose)) {
                message("Returning CPM from pre-calculated 'cpm' assay.")
            }
            return(assay(x = object, i = "cpm"))
        }
        ## Otherwise, calculate on the fly.
        assert(is.numeric(sizeFactors(object)))
        if (isTRUE(verbose)) {
            message("Calculating CPM with 'scater::calculateCPM()'.")
        }
        calculateCPM(object)
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
    function(object, assay = NULL, verbose = FALSE) {
        assert(isFlag(verbose))
        ## Check for pre-calculated CPM (not typical).
        method <- .seuratNormalizationMethod(object, assay = assay)
        scaleFactor <- .seuratScaleFactor(object, assay = assay)
        if (!(method == "RC" && scaleFactor == 1e6L)) {
            if (isTRUE(verbose)) {
                message("Generating CPM with 'Seurat::NormalizeData()'.")
            }
            object <- NormalizeData(
                object = object,
                assay = assay,
                normalization.method = "RC",
                scale.factor = 1e6L,
                verbose = verbose
            )
        }
        if (isTRUE(verbose)) {
            message("Returning CPM with 'Seurat::GetAssayData()'.")
        }
        GetAssayData(object = object, slot = "data", assay = assay)
    }



#' @rdname cpm
#' @export
setMethod(
    f = "cpm",
    signature = signature("Seurat"),
    definition = `cpm,Seurat`
)
