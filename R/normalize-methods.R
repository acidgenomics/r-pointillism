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
#' @note Updated 2019-10-26.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @return Modified object.
#'
#' @seealso
#' - [estimateSizeFactors()].
#' - [SingleCellExperiment::normcounts()].
#' - [SingleCellExperiment::logcounts()].
#' - [scater::normalizeSCE()].
#' - [Seurat::NormalizeData()].
#' - [monocle3::preprocess_cds()].
#' - [monocle3::normalized_counts()].
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
#' object <- normalize(object)
#'
#' ## Seurat ====
#' object <- Seurat
#' object <- normalize(object)
#'
#' ## cell_data_set ====
#' ## > object <- cell_data_set
#' ## > object <- normalize(object)
NULL



#' @rdname normalize
#' @name normalize
#' @importFrom BiocGenerics normalize
#' @usage normalize(object, ...)
#' @export
NULL



## Updated 2019-10-30.
`normalize,SingleCellExperiment` <-  # nolint
    function(object, verbose = FALSE) {
        validObject(object)
        assert(isFlag(verbose))
        if (is.null(sizeFactors(object))) {
            object <- estimateSizeFactors(object)
        }
        assert(!is.null(sizeFactors(object)))
        if (isTRUE(verbose)) {
            message(
                "Computing 'normcounts' and 'logcounts' assays using ",
                "'scater::normalizeSCE()'."
            )
        }
        ## Get normcounts assay.
        sce <- normalizeSCE(
            object = object,
            return_log = FALSE,
            ## centreSizeFactors is now deprecated.
            centre_size_factors = FALSE
        )
        assert(
            isSubset("normcounts", assayNames(sce)),
            !isSubset("logcounts", assayNames(sce))
        )
        normcounts <- normcounts(sce)
        ## Get logcounts assay.
        sce <- normalizeSCE(
            object = object,
            return_log = TRUE,
            centre_size_factors = FALSE
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



## Updated 2019-08-05.
`normalize,Seurat` <-  # nolint
    function(object, verbose = FALSE) {
        assert(isFlag(verbose))
        if (isTRUE(verbose)) {
            message("Normalizing with 'Seurat::NormalizeData()'.")
        }
        NormalizeData(object = object, verbose = verbose)
    }



#' @rdname normalize
#' @export
setMethod(
    f = "normalize",
    signature = signature("Seurat"),
    definition = `normalize,Seurat`
)



## Updated 2019-08-05.
`normalize,cell_data_set` <-  # nolint
    function(object, verbose = FALSE) {
        assert(isFlag(verbose))
        if (isTRUE(verbose)) {
            message("Normalizing with 'monocle3::preprocess_cds()'.")
        }
        monocle3::preprocess_cds(
            cds = object,
            method = "PCA",
            norm_method = "log",
            pseudo_count = 1L,
            scaling = TRUE,
            verbose = verbose
        )
    }



#' @rdname normalize
#' @export
setMethod(
    f = "normalize",
    signature = signature("cell_data_set"),
    definition = `normalize,cell_data_set`
)
