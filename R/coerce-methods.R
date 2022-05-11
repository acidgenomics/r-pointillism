#' Force an object to belong to a class
#'
#' @name coerce
#' @importFrom methods coerce
#' @exportMethod coerce
#' @note Updated 2022-05-11.
#'
#' @examples
#' data(seurat)
#' data(SingleCellExperiment_Seurat, package = "AcidTest")
#'
#' ## SingleCellExperiment to Seurat ====
#' from <- SingleCellExperiment_Seurat
#' to <- as(from, "Seurat")
#' class(to)
#' print(to)
#'
#' ## Seurat to SingleCellExperiment ====
#' from <- seurat
#' to <- as(from, "SingleCellExperiment")
#' class(to)
#' print(to)
NULL



## Updated 2022-05-11.
`coerce,Seurat,SCE` <-  # nolint
    function(from) {
        as.SingleCellExperiment(from)
    }

## Updated 2022-05-11.
`coerce,Seurat,RSE` <-  # nolint
    function(from) {
        to <- from
        to <- as(to, "SingleCellExperiment")
        to <- as(to, "RangedSummarizedExperiment")
        to
    }

## Updated 2022-05-11.
`coerce,Seurat,SE` <-  # nolint
    function(from) {
        to <- from
        to <- as(to, "RangedSummarizedExperiment")
        to <- as(to, "SummarizedExperiment")
        to
    }

## Updated 2022-05-11.
`coerce,SCE,Seurat` <-  # nolint
    function(from) {
        as.Seurat(from)
    }



#' @rdname coerce
#' @name coerce,Seurat,SingleCellExperiment-method
setAs(
    from = "Seurat",
    to = "SingleCellExperiment",
    def = `coerce,Seurat,SCE`
)

#' @rdname coerce
#' @name coerce,Seurat,RangedSummarizedExperiment-method
setAs(
    from = "Seurat",
    to = "RangedSummarizedExperiment",
    def = `coerce,Seurat,RSE`
)

#' @rdname coerce
#' @name coerce,Seurat,SummarizedExperiment-method
setAs(
    from = "Seurat",
    to = "SummarizedExperiment",
    def = `coerce,Seurat,SE`
)

#' @rdname coerce
#' @name coerce,SingleCellExperiment,Seurat-method
setAs(
    from = "SingleCellExperiment",
    to = "Seurat",
    def = `coerce,SCE,Seurat`
)
