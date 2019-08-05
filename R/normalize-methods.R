## Improve spike-in control handling support here, similar to scater.

## FIXME `log_exprs_offset`



#' Normalize expression using pre-computed size factors
#'
#' This function calculates size factor normalized and log normalized counts
#' from the raw count matrix defined in [counts()]. These matrices are stored in
#' [normcounts()] and [logcounts()] slots in [assays()], respectively.
#'
#' If no library size factors are defined in [sizeFactors()], these values will
#' be computed internally automatically using [estimateSizeFactors()] with the
#' recommended default settings.
#'
#' For complex normalizations involving custom size factors or spike-ins (i.e.
#' when [spikeNames()]) is defined, call [scater::normalizeSCE()] directly
#' instead.
#'
#' @section Normalized counts:
#'
#' Normalized counts are computed by dividing the counts for each cell by the
#' size factor for that cell. This aims to remove cell-specific scaling biases,
#' due to differences in sequencing coverage or capture efficiency.
#'
#' @section Log normalized counts:
#'
#' Log-normalized values are calculated by adding a pseudocount offset of 1 to
#' the normalized count and performing a [`log2()`][base::log2] transformation.
#'
#' @section Centering at unity:
#'
#' When centering is applied to the size factors (recommended by default), all
#' sets of size factors will be adjusted to have the same mean prior to
#' calculation of normalized expression values. This ensures that abundances are
#' roughly comparable between features normalized with different sets of size
#' factors. By default, the center mean is unity, which means that the computed
#' expression values can be interpreted as being on the same scale as the
#' log2-counts. It also means that the added offset (i.e. `1`) added to the
#' normalized counts during the log transformation step can be interpreted as a
#' pseudo-count (i.e., on the same scale as the counts).
#'
#' @name normalize
#' @note Updated 2019-08-05.
#'
#' @seealso
#' - `estimateSizeFactors()`.
#' - `SingleCellExperiment::normcounts()`.
#' - `SingleCellExperiment::logcounts()`.
#' - `scater::normalizeSCE()`.
#' - `Seurat::NormalizeData()`.
#' - `monocle3::normalized_counts()`.
#'
#' @examples
#' data(SingleCellExperiment, package = "acidtest")
#'
#' ## SingleCellExperiment ====
#' object <- SingleCellExperiment
#' object <- normalize(object)
#' class(normcounts(object))
#' class(logcounts(object))
NULL



## Updated 2019-08-05.
`normalize,SingleCellExperiment` <-  # nolint
    function(object) {
        validObject(object)
        ## Use `scater::normalizeSCE()` instead when working with spike-ins.
        assert(!hasLength(spikeNames(object)))
        if (is.null(sizeFactors(object))) {
            object <- estimateSizeFactors(object)
        }
        assert(!is.null(sizeFactors(object)))
        message(paste(
            "Computing normcounts and logcounts using",
            "`scater::normalizeSCE()`."
        ))

        ## Shared arguments for `normalizeSCE()` calls.
        args <- list(
            object = object,
            exprs_values = "counts",
            centre_size_factors = TRUE,
            preserve_zeroes = FALSE
        )

        ## Get normcounts assay.
        sce <- do.call(
            what = normalizeSCE,
            args = c(
                args,
                return_log = FALSE
            )
        )
        assert(
            isSubset("normcounts", assayNames(sce)),
            !isSubset("logcounts", assayNames(sce))
        )
        normcounts <- normcounts(sce)

        ## Get logcounts assay.
        sce <- do.call(
            what = normalizeSCE,
            args = c(
                args,
                return_log = TRUE,
                log_exprs_offset = 1L
            )
        )
        assert(
            isSubset("logcounts", assayNames(sce)),
            !isSubset("normcounts", assayNames(sce))
        )
        logcounts <- logcounts(sce)

        ## Slot the normalized counts in object.
        normcounts(object) <- normcounts
        logcounts(object) <- logcounts

        ## Stash scater package version in metadata.
        metadata(object)[["scater"]] <- packageVersion("scater")

        object
    }



#' @rdname normalize
#' @export
setMethod(
    f = "normalize",
    signature = signature("SingleCellExperiment"),
    definition = `normalize,SingleCellExperiment`
)
