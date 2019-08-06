#' Size factor adjusted normalized counts
#'
#' @name normcounts
#' @keywords internal
#' @note Updated 2019-08-05.
#'
#' @inheritParams acidroxygen::params
#'
#' @return `sparseMatrix`.
#'
#' @seealso
#' - [normalize()].
#' - `scater::normalizeSCE()`.
#' - `scater::normalizeCounts()`.
#' - `Seurat::NormalizeData()`.
#' - `monocle3::prepare_cds()`.
#' - `monocle3::normalized_counts()`.
#'
#' @examples
#' data(
#'     Seurat,
#'     SingleCellExperiment,
#'     cell_data_set,
#'     package = "acidtest"
#' )
#'
#' ## SingleCellExperiment ===
#' object <- SingleCellExperiment
#' object <- normalize(object)
#' normcounts <- normcounts(object)
#' class(normcounts)
#'
#' ## Seurat ====
#' object <- Seurat
#' normcounts <- normcounts(object)
#' class(normcounts)
#'
#' ## cell_data_set ====
#' object <- cell_data_set
#' object <- normalize(object)
#' normcounts <- normcounts(object)
#' class(normcounts)
NULL



#' @rdname normcounts
#' @name normcounts
#' @importFrom SingleCellExperiment normcounts
#' @usage normcounts(object, ...)
#' @export
NULL



## Updated 2019-08-05.
`normcounts,Seurat` <-  # nolint
    function(object, assay = NULL, verbose = FALSE) {
        assert(isFlag(verbose))
        ## Check for pre-calculated relative counts (not typical).
        method <- .seuratNormalizationMethod(object, assay = assay)
        scaleFactor <- .seuratScaleFactor(object, assay = assay)
        if (!(method == "RC" && scaleFactor == 1L)) {
            if (isTRUE(verbose)) {
                message(paste(
                    "Generating normalized counts with",
                    "`Seurat::NormalizeData()`."
                ))
            }
            object <- NormalizeData(
                object = object,
                normalization.method = "RC",
                scale.factor = 1L,
                verbose = verbose
            )
        }
        if (isTRUE(verbose)) {
            message(paste(
                "Returning normalized counts with",
                "`Seurat::GetAssayData()`."
            ))
        }
        GetAssayData(object = object, assay = assay, slot = "data")
    }



#' @rdname normcounts
#' @export
setMethod(
    f = "normcounts",
    signature = signature("Seurat"),
    definition = `normcounts,Seurat`
)



## Updated 2019-08-05.
`normcounts,cell_data_set` <-  # nolint
    function(object, verbose = FALSE) {
        assert(isFlag(verbose))
        if (isTRUE(verbose)) {
            message(paste(
                "Getting normalized counts with",
                "`monocle3::normalized_counts()`."
            ))
        }
        monocle3::normalized_counts(
            cds = object,
            norm_method = "size_only",
            pseudocount = NULL
        )
    }



#' @rdname normcounts
#' @export
setMethod(
    f = "normcounts",
    signature = signature("cell_data_set"),
    definition = `normcounts,cell_data_set`
)
