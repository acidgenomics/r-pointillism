#' Log normalized counts
#'
#' @name logcounts
#' @note Updated 2019-08-05.
#'
#' @inheritParams acidroxygen::params
#'
#' @return `sparseMatrix`.
#'
#' @seealso
#' - [normalize()].
#' - `scater::normalizeSCE()`.
#' - `Seurat::NormalizeData()`.
#' - `monocle3::prepare_cds()`.
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
#' logcounts <- logcounts(object)
#' class(logcounts)
#'
#' ## Seurat ====
#' object <- Seurat
#' logcounts <- logcounts(object)
#' class(logcounts)
#'
#' ## cell_data_set ====
#' object <- cell_data_set
#' object <- normalize(object)
#' logcounts <- logcounts(object)
#' class(logcounts)
NULL



#' @rdname logcounts
#' @name logcounts
#' @importFrom SingleCellExperiment logcounts
#' @usage logcounts(object, ...)
#' @export
NULL



## Updated 2019-08-05.
`logcounts,Seurat` <-  # nolint
    function(object, assay = NULL, verbose = FALSE) {
        assert(isFlag(verbose))
        norm <- .seuratNormalizationMethod(object, assay = assay)
        if (norm != "LogNormalize") {
            if (isTRUE(verbose)) {
                message(paste(
                    "Generating log normalized counts with",
                    "`Seurat::NormalizeData()`."
                ))
            }
            object <- NormalizeData(
                object = object,
                normalization.method = "LogNormalize",
                verbose = verbose
            )
        }
        if (isTRUE(verbose)) {
            message(paste(
                "Returning log normalized counts with",
                "`Seurat::GetAssayData()`."
            ))
        }
        GetAssayData(object = object, assay = assay, slot = "data")
    }



#' @rdname counts
#' @export
setMethod(
    f = "logcounts",
    signature = signature("Seurat"),
    definition = `logcounts,Seurat`
)



## Updated 2019-08-05.
`logcounts,cell_data_set` <-  # nolint
    function(object, verbose = FALSE) {
        assert(isFlag(verbose))
        if (isTRUE(verbose)) {
            message(paste(
                "Returning log normalized counts with",
                "`monocle3::normalized_counts()`."
            ))
        }
        monocle3::normalized_counts(
            cds = object,
            norm_method = "log",
            pseudocount = 1L
        )
    }



#' @rdname counts
#' @export
setMethod(
    f = "logcounts",
    signature = signature("cell_data_set"),
    definition = `logcounts,cell_data_set`
)
