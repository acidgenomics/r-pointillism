#' Force an object to belong to a class
#'
#' @name coerce
#' @importFrom methods coerce
#' @exportMethod coerce
#' @note Updated 2019-07-31.
#'
#' @seealso
#' - `Seurat::CreateSeuratObject()`.
#' - `Seurat::as.Seurat()`.
#' - `Seurat::as.SingleCellExperiment()`.
#'
#' @examples
#' data(Seurat, SingleCellExperiment, package = "acidtest")
#'
#' ## SingleCellExperiment to Seurat ====
#' x <- as(SingleCellExperiment, "Seurat")
#' class(x)
#' print(x)
#'
#' ## Seurat to SingleCellExperiment ====
#' x <- as(Seurat, "SingleCellExperiment")
#' print(x)
NULL



## Updated 2019-07-31.
`coerce,Seurat,SingleCellExperiment` <-  # nolint
    function(from) {
        ## Using the Seurat S3 coercion method here.
        to <- as.SingleCellExperiment(x = from, assay = NULL)
        ## Row and column data.
        rowRanges(to) <- rowRanges(from)
        colData(to) <- colData(to)
        ## Metadata.
        metadata(to) <- metadata(from)
        metadata(to)[["scaleData"]] <- GetAssayData(from, slot = "scale.data")
        metadata(to)[["variableFeatures"]] <- VariableFeatures(from)
        to
    }



#' @rdname coerce
#' @name coerce,Seurat,SingleCellExperiment-method
#'
#' @section `Seurat` to `SingleCellExperiment`:
#' S4 coercion support for creating a `SingleCellExperiment` from a `Seurat`
#' class object. The [Seurat FAQ page](https://satijalab.org/Seurat/faq)
#' explains the `Seurat` S4 class structure in detail. Internally, this method
#' improves the basic `Seurat::as.SingleCellExperiment` S3 coercion method,
#' including the `object@scale.data` matrix, and will keep track of stashed
#' `rowRanges` and `metadata` if the `Seurat` object was originally created
#' from a `SingleCellExperiment` (i.e. from the bcbioSingleCell package).
setAs(
    from = "Seurat",
    to = "SingleCellExperiment",
    def = `coerce,Seurat,SingleCellExperiment`
)



## Updated 2019-07-31.
`coerce,Seurat,RangedSummarizedExperiment` <-  # nolint
    function(from) {
        to <- from
        to <- as(to, "SingleCellExperiment")
        to <- as(to, "RangedSummarizedExperiment")
        to
    }



#' @rdname coerce
#' @name coerce,Seurat,RangedSummarizedExperiment-method
#'
#' @section `Seurat` to `RangedSummarizedExperiment`:
#' S4 coercion support for creating a `RangedSummarizedExperiment` from a
#' `Seurat` class object.
setAs(
    from = "Seurat",
    to = "RangedSummarizedExperiment",
    def = `coerce,Seurat,RangedSummarizedExperiment`
)



## Updated 2019-07-31.
`coerce,Seurat,SummarizedExperiment` <-  # nolint
    function(from) {
        to <- from
        to <- as(to, "RangedSummarizedExperiment")
        to <- as(to, "SummarizedExperiment")
        to
    }



#' @rdname coerce
#' @name coerce,Seurat,SummarizedExperiment-method
#' @section `Seurat` to `SummarizedExperiment`:
#' S4 coercion support for creating a `SummarizedExperiment` from a `Seurat`
#' class object.
setAs(
    from = "Seurat",
    to = "SummarizedExperiment",
    def = `coerce,Seurat,SummarizedExperiment`
)



## Updated 2019-07-31.
`coerce,SingleCellExperiment,Seurat` <-  # nolint
    function(from) {
        ## Create the Seurat object. Note that `as.Seurat()` method requires
        ## `logcounts` to be defined in `assays()`, so we're using
        ## `CreateSeuratObject()` here instead.
        to <- CreateSeuratObject(
            counts = counts(from),
            project = "pointillism",
            assay = "RNA",
            ## Already applied filtering cutoffs for cells and genes.
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
        metadata[["sessionInfo"]] <- session_info()
        ## Seurat v3 still recommends using `misc` slot.
        misc <- list(
            rowRanges = rowRanges,
            metadata = metadata
        )
        misc <- Filter(Negate(is.null), misc)
        slot(to, name = "misc") <- misc
        ## Return.
        to
    }



#' @rdname coerce
#' @name coerce,SingleCellExperiment,Seurat-method
#'
#' @section `SingleCellExperiment` to `Seurat`:
#' Interally `Seurat::CreateSeuratObject` is called without applying any
#' additional filtering cutoffs, since we have already defined them during our
#' quality control analysis. Here we are passing the raw gene-level counts of
#' the filtered cells into a new `Seurat` class object. Use
#' `convertGenesToSymbols` to convert gene IDs to names (symbols).
setAs(
    from = "SingleCellExperiment",
    to = "Seurat",
    def = `coerce,SingleCellExperiment,Seurat`
)
