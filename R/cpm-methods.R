## FIXME Consider moving this to AcidSingleCell.
## FIXME If that is the case, scater needs to be added as a suggested package.



#' @name cpm
#' @inherit AcidGenerics::cpm
#' @keywords internal
#' @note Updated 2020-01-30.
#'
#' @inheritParams AcidRoxygen::params
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
#'     package = "AcidTest"
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
NULL



## Updated 2020-01-30.
`cpm,SCE` <-  # nolint
    function(object) {
        ## Early return if cpm assay is defined.
        if (isSubset("cpm", assayNames(object))) {
            return(assay(x = object, i = "cpm"))
        }
        ## Otherwise, calculate on the fly.
        assert(is.numeric(sizeFactors(object)))
        alert("Calculating CPM with {.pkg scater}::{.fun calculateCPM}.")
        calculateCPM(object)
    }



## Updated 2020-01-30.
`cpm,Seurat` <-  # nolint
    function(object, assay = NULL) {
        ## Check for pre-calculated CPM (not typical).
        method <- .seuratNormalizationMethod(object, assay = assay)
        scaleFactor <- .seuratScaleFactor(object, assay = assay)
        if (!(method == "RC" && scaleFactor == 1e6L)) {
            alert(
                "Generating CPM with {.pkg Seurat}::{.fun NormalizeData}."
            )
            object <- NormalizeData(
                object = object,
                assay = assay,
                normalization.method = "RC",
                scale.factor = 1e6L,
                verbose = TRUE
            )
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

#' @rdname cpm
#' @export
setMethod(
    f = "cpm",
    signature = signature("SingleCellExperiment"),
    definition = `cpm,SCE`
)
