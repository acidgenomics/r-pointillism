#' Estimate size factors
#'
#' Define size factors from the library sizes, and then apply centering at
#' unity. This ensures that the library size adjustment yields values comparable
#' to those generated after normalization with other sets of size factors.
#'
#' The estimated size factors computed by this function can be accessed using
#' the accessor function [sizeFactors()]. Alternatively library size estimators
#' can also be supplied using the assignment function [sizeFactors<-()].
#'
## Note that we're computing internally on the count matrix as a DelayedArray,
## so we can handle millions of cells without the calculations blowing up in
## memory.
#'
#' @name estimateSizeFactors
#' @note Updated 2019-08-05.
#'
#' @inheritParams params
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
#' @param center `numeric(1)`.
#'   If non-zero, scales all size factors so that the average size factor across
#'   cells is equal to the value defined. Set to `0` to disable centering.
#'
#' @param ... Additional arguments.
#'
#' @seealso
#' From DESeq2:
#' - `DESeq2::estimateSizeFactors()`.
#' - `DESeq2::estimateSizeFactorsForMatrix().`
#'
#' From scater:
#' - `scater::librarySizeFactors()`.
#' - `scater::centreSizeFactors()`.
#' - `scater::normalizeSCE()`.
#'
#' From monocle3:
#' - `monocle3::estimate_size_factors()`.
#' - `monocle3:::estimate_sf_sparse()`.
#'
#' @return Modified object.
#'   Use `[sizeFactors()] to access the computed size factor numeric.
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



## Updated 2019-08-05.
.librarySizeFactors <-  # nolint
    function(
        counts,
        type = c(
            "mean-ratio",
            "geometric-mean-ratio",
            "log-geometric-mean-ratio",
            "deseq2-median-ratio"
        )
    ) {
        assert(
            is(counts, "DelayedArray"),
            !anyNA(counts)
        )
        type <- match.arg(type)
        message(sprintf(
            fmt = "Calculating library size factors using %s method.",
            type
        ))

        ## Get the sum of expression per cell.
        libSizes <- colSums2(counts)

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
                DESeq2::estimateSizeFactorsForMatrix(counts)
            }
        )

        assert(!anyNA(sf))
        names(sf) <- colnames(counts)

        sf
    }



## Updated 2019-08-05.
.centerSizeFactors <- function(sf, center = 1L) {
    assert(
        is.numeric(sf),
        hasNames(sf),
        isScalarNumeric(center),
        isPositive(center)
    )
    message(sprintf("Centering size factors at %d.", center))
    sf <- sf / mean(sf) * center
    assert(mean(sf) == center)
    sf
}



## Updated 2019-08-05.
`estimateSizeFactors,SingleCellExperiment` <-  # nolint
    function(object, type, center
    ) {
        validObject(object)
        assert(
            isScalarNumeric(center),
            isNonNegative(center)
        )
        type <- match.arg(type)
        counts <- DelayedArray(counts(object))
        assert(is(counts, "DelayedMatrix"))
        sf <- .librarySizeFactors(counts = counts, type = type)
        if (center > 0L) {
            sf <- .centerSizeFactors(sf = sf, center = center)
        }
        assert(is.numeric(sf), hasNames(sf))
        sizeFactors(object, type = NULL) <- sf
        object
    }

formals(`estimateSizeFactors,SingleCellExperiment`)[
    c("type", "center")
] <- list(
    type = formals(.librarySizeFactors)[["type"]],
    center = formals(.centerSizeFactors)[["center"]]
)



#' @rdname estimateSizeFactors
#' @export
setMethod(
    f = "estimateSizeFactors",
    signature = signature("SingleCellExperiment"),
    definition = `estimateSizeFactors,SingleCellExperiment`
)



## Updated 2019-08-05.
`estimateSizeFactors,cell_data_set` <-  # nolint
    `estimateSizeFactors,SingleCellExperiment`



#' @rdname estimateSizeFactors
#' @export
setMethod(
    f = "estimateSizeFactors",
    signature = signature("cell_data_set"),
    definition = `estimateSizeFactors,cell_data_set`
)
