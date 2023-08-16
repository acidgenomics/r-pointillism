#' @name as.SingleCellExperiment
#' @inherit AcidGenerics::as.SingleCellExperiment
#' @note Updated 2023-08-16.
#'
#' @details
#' S4 coercion support for creating a `SingleCellExperiment` from a `Seurat`
#' class object. Internally, this method improves the basic
#' `Seurat::as.SingleCellExperiment` S3 coercion method, including the
#' `object@scale.data` matrix, and will keep track of stashed `rowRanges` and
#' `metadata` if the `Seurat` object was originally created from a
#' `SingleCellExperiment` (i.e. from the bcbioSingleCell package).
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @seealso
#' - `Seurat::as.SingleCellExperiment()`.
#'
#' @examples
#' data(Seurat, package = "AcidTest")
#'
#' ## Seurat to SingleCellExperiment ====
#' from <- Seurat
#' to <- as.SingleCellExperiment(from)
#' class(to)
#' print(to)
NULL



## FIXME Seems like as.SummarizedExperiment method in Seurat 5 is dropping
## scale.data.

## Updated 2023-08-16.
`as.SingleCellExperiment,Seurat` <- # nolint
    function(x) {
        from <- x
        validObject(from)
        to <- Seurat::as.SingleCellExperiment(x = from, assay = NULL)
        if (isSubset("ident", colnames(colData(to)))) {
            idents <- Idents(from)
            assert(
                identical(names(idents), rownames(colData(to))),
                is.factor(idents),
                !all(is.na(idents))
            )
            colData(to)[["ident"]] <- unname(idents)
        }
        rowRanges(to) <- rowRanges(from)
        metadata(to) <- metadata(from)
        ## FIXME scale.data is no longer defined in Seurat 5?
        ## FIXME Only do this conditionally if defined.
        metadata(to)[["scaleData"]] <- LayerData(from, layer = "scale.data")
        metadata(to)[["variableFeatures"]] <- VariableFeatures(from)
        validObject(to)
        to
    }



#' @rdname as.SingleCellExperiment
#' @export
setMethod(
    f = "as.SingleCellExperiment",
    signature = signature(
        x = "Seurat"
    ),
    definition = `as.SingleCellExperiment,Seurat`
)
