## Monocle guide mentions using DelayedArray approach for this.



#' Estimate size factors
#'
#' @name estimateSizeFactors
#' @inheritParams params
#'
#' @seealso
#' - `DESeq2::estimateSizeFactors()`.
#' - `monocle3::estimate_size_factors()`.
#'
#' @examples
#' data(
#'     SingleCellExperiment,
#'     cell_data_set,
#'     package = "acidtest"
#' )
#'
#' ## SingleCellExperiment ====
#' object <- SingleCellExperiment
#' sizeFactors(object)
#' object <- estimateSizeFactors(object)
#' sizeFactors(object)
#'
#' ## cell_data_set ====
#' object <- cell_data_set
#' object <- estimateSizeFactors(cell_data_set)
#' sizeFactors(object)
NULL



## Consider porting some of the DESeqDataSet method code here.
## Updated 2019-08-02.
`estimateSizeFactors,SingleCellExperiment` <-  # nolint
    function(object) {
        colData(object)[["sizeFactor"]] <- NULL
        cds <- as(object, "cell_data_set")
        cds <- estimateSizeFactors(cds)
        assert(
            is(cds, "cell_data_set"),
            identical(dimnames(object), dimnames(cds)),
            isSubset("Size_Factor", colnames(colData(cds)))
        )
        sf <- colData(cds)[["Size_Factor"]]
        assert(is.numeric(sf))
        colData(object)[["sizeFactor"]] <- unname(sf)
        object
    }



#' @rdname estimateSizeFactors
#' @export
setMethod(
    f = "estimateSizeFactors",
    signature = signature("SingleCellExperiment"),
    definition = `estimateSizeFactors,SingleCellExperiment`
)



## This method sparse calculations, which is awesome.
## Updated 2019-08-02.
`estimateSizeFactors,cell_data_set` <-  # nolint
    function(object) {
        monocle3::estimate_size_factors(
            cds = object,
            round_exprs = TRUE,
            method = "mean-geometric-mean-total"
        )
    }



#' @rdname estimateSizeFactors
#' @export
setMethod(
    f = "estimateSizeFactors",
    signature = signature("cell_data_set"),
    definition = `estimateSizeFactors,cell_data_set`
)
