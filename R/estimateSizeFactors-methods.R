## Monocle3 mentions using DelayedArray approach for this. Consider adapting
## code to use this package, which is capable of handling millions of cells.

## FIXME Refer to DESeqDataSet approach here for ideas.
## FIXME Consider not exporting colData for Seurat...use internally?
## FIXME Cover that this returns identical values as monocle3 function.



#' Estimate size factors
#'
#' @name estimateSizeFactors
#' @note Updated 2019-08-04.
#'
#' @inheritParams params
#' @param colname `character(1)`.
#'   Column name where to define size factor values in
#'   [`colData()`][SummarizedExperiment::colData].
#' @param roundExprs `logical(1)`.
#'   Round expression values to integer counts, prior to calculation.
#' @param method `character(1)`.
#'   Size factor normalization method.
#' @param ... Additional arguments.
#'
#' @seealso
#' - `DESeq2::estimateSizeFactors()`.
#' - `monocle3::estimate_size_factors()`.
#'
#' @return
#' - `Matrix`: Named `numeric`, containing size factor values.
#' - `SingleCellExperiment`: Modified object, containing size factor values
#'   as `sizeFactor` column in [`colData()`][SummarizedExperiment::colData].
#' - `cell_data_set`: Modified object, containing size factor values
#'   as `Size_Factor` column in [`colData()`][SummarizedExperiment::colData].
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



#' @rdname estimateSizeFactors
#' @name estimateSizeFactors
#' @importFrom BiocGenerics estimateSizeFactors
#' @usage estimateSizeFactors(object, ...)
#' @export
NULL



## Modified version from monocle3:
## - `monocle3::estimate_size_factors()`.
## - `monocle3:::estimate_sf_dense()`.
## - `monocle3:::estimate_sf_sparse()`.
## Updated 2019-08-04.
`estimateSizeFactors,Matrix`  <-  # nolint
    function(
        object,
        roundExprs = FALSE,
        method = c(
            "mean-geometric-mean-total",
            "mean-geometric-mean-log-total"
        )
    ) {
        assert(
            !anyNA(object),
            isFlag(roundExprs)
        )
        method <- match.arg(method)
        message(sprintf(
            fmt = "Calculating size factors using %s method.",
            method
        ))

        ## Round the counts to integer values, if desired.
        if (isTRUE(roundExprs)) {
            object <- round(object, digits = 0L)
        }

        ## Get the sum of expression per cell.
        sum <- Matrix::colSums(object)

        ## Check for zero count cells and inform the user where in matrix.
        zero <- sum == 0L
        if (isTRUE(any(zero))) {
            stop(sprintf(
                fmt = "Cells with no expression detected: %s",
                toString(
                    unname(which(zero)),
                    width = 100L
                )
            ))
        }

        ## Calculate the size factors per cell.
        sf <- switch(
            EXPR = method,
            "mean-geometric-mean-total" = {
                sum / exp(mean(log(sum)))
            },
            "mean-geometric-mean-log-total" = {
                log(sum) / exp(mean(log(log(sum))))
            }
        )
        sf[is.na(sf)] <- 1L
        sf
    }



#' @rdname estimateSizeFactors
#' @export
setMethod(
    f = "estimateSizeFactors",
    signature = signature("Matrix"),
    definition = `estimateSizeFactors,Matrix`
)



## Consider porting some of the DESeqDataSet method code here.
## Updated 2019-08-04.
`estimateSizeFactors,SingleCellExperiment` <-  # nolint
    function(object, colname = "sizeFactor") {
        validObject(object)
        assert(isString(colname))
        if (isSubset(colname, colnames(colData(object)))) {
            message(sprintf(
                fmt = paste(
                    "Replacing size factor values defined in",
                    "`%s` column of `colData()`."
                ),
                colname
            ))
        }
        counts <- counts(object)
        sf <- estimateSizeFactors(counts)
        colData(object)[[colname]] <- unname(sf)
        object
    }



#' @rdname estimateSizeFactors
#' @export
setMethod(
    f = "estimateSizeFactors",
    signature = signature("SingleCellExperiment"),
    definition = `estimateSizeFactors,SingleCellExperiment`
)



## Updated 2019-08-04.
`estimateSizeFactors,cell_data_set` <-  # nolint
    `estimateSizeFactors,SingleCellExperiment`

formals(`estimateSizeFactors,cell_data_set`)[["colname"]] <- "Size_Factor"



#' @rdname estimateSizeFactors
#' @export
setMethod(
    f = "estimateSizeFactors",
    signature = signature("cell_data_set"),
    definition = `estimateSizeFactors,cell_data_set`
)
