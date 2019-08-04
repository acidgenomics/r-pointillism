#' Estimate size factors
#'
#' @name estimateSizeFactors
#' @note Updated 2019-08-04.
#'
#' @inheritParams params
#' @param colname `character(1)`.
#'   Column name where to define size factor values in
#'   [`colData()`][SummarizedExperiment::colData].
#' @param type `character(1)`.
#'   Method for estimation:
#'   ```
#'   libSize <- colSums(counts(object))
#'   ```
#'   - `mean-ratio`: From scater.
#'     ```
#'     libSize / mean(libSize)
#'     ```
#'   - `geometric-mean-ratio`: From monocle3.
#'     ```
#'     libSize / geometricMean(libSize)
#'     ```
#'   - `mean-geometric-mean-log-total`: From monocle3.
#'     ```
#'     log(libSize) / geometricMean(log(libSize))
#'     ```
#'   - `median-ratio`: From DESeq2.
#'     The size factor is a median ratio of the sample over a "pseudosample":
#'     for each feature (i.e. gene), the geometric mean of all samples.
#'
#' @param ... Additional arguments.
#'
#' @seealso
#' - `scater::librarySizeFactors()`.
#' - `monocle3::estimate_size_factors()`.
#' - `monocle3:::estimate_sf_sparse()`.
#' - `DESeq2::estimateSizeFactors()`.
#' - `DESeq2::estimateSizeFactorsForMatrix().`
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



## Updated 2019-08-04.
`estimateSizeFactors,DelayedArray`  <-  # nolint
    function(
        object,
        type = c(
            "mean-ratio",
            "geometric-mean-ratio",
            "log-geometric-mean-ratio",
            "deseq2-median-ratio"
        )
    ) {
        assert(!anyNA(object))
        type <- match.arg(type)
        message(sprintf("Calculating size factors using %s method.", type))

        ## Get the sum of expression per cell.
        libSizes <- colSums2(object)

        ## Error on detection of cells without any expression.
        zero <- libSizes == 0L
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
            EXPR = type,
            "mean-ratio" = {
                libSizes / mean(libSizes)
            },
            "geometric-mean-ratio" = {
                libSizes / geometricMean(libSizes)
            },
            "log-geometric-mean-ratio" = {
                log(libSizes) / geometricMean(log(libSizes))
            },
            "deseq2-median-ratio" = {
                estimateSizeFactorsForMatrix(object)
            }
        )

        assert(!anyNA(sf))
        names(sf) <- colnames(object)

        sf
    }



#' @rdname estimateSizeFactors
#' @export
setMethod(
    f = "estimateSizeFactors",
    signature = signature("DelayedArray"),
    definition = `estimateSizeFactors,DelayedArray`
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
        counts <- DelayedArray(counts(object))
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
