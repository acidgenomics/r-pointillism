#' @name as.Seurat
#' @inherit AcidGenerics::as.Seurat
#' @note Updated 2022-05-11.
#'
#' @details
#' Note that `Seurat::as.Seurat()` method requires `logcounts` to be defined
#' in `assays()`, so we're using `CreateSeuratObject()` here instead.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @seealso
#' - `Seurat::CreateSeuratObject()`.
#' - `Seurat::as.Seurat()`.
#'
#' @examples
#' data(SingleCellExperiment_Seurat, package = "AcidTest")
#'
#' ## SingleCellExperiment to Seurat ====
#' from <- SingleCellExperiment_Seurat
#' to <- as.Seurat(from)
#' class(to)
#' print(to)
NULL



#' Coerce SCE to Seurat
#'
#' @note Updated 2022-05-11.
#' @noRd
#'
#' @details
#' Note that `Seurat::as.Seurat()` method requires `logcounts` to be defined
#' in `assays()`, so we're using `CreateSeuratObject()` here instead.
`as.Seurat,SCE` <-  # nolint
    function(x) {
        from <- x
        to <- CreateSeuratObject(
            counts = counts(from),
            project = "pointillism",
            assay = "RNA",
            min.cells = 0L,
            min.features = 0L,
            names.field = 1L,
            names.delim = "_",
            meta.data = as.data.frame(colData(from))
        )
        ## Check that the dimensions match exactly.
        assert(identical(x = dim(from), y = dim(to)))
        ## Stash rowRanges.
        rowRanges <- rowRanges(from)
        ## Stash metadata.
        metadata <- metadata(from)
        ## Update the session information.
        metadata[["sessionInfo"]] <- sessionInfo()
        ## Seurat v3 still recommends using `misc` slot.
        misc <- list(
            "rowRanges" = rowRanges,
            "metadata" = metadata
        )
        misc <- Filter(Negate(is.null), misc)
        slot(to, name = "misc") <- misc
        to
    }



#' @rdname as.Seurat
#' @export
setMethod(
    f = "as.Seurat",
    signature = signature(
        x = "SingleCellExperiment"
    ),
    definition = `as.Seurat,SCE`
)
