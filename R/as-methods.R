#' @inherit methods::as
#' @name as
#' @aliases coerce
#' @importFrom methods coerce
#' @exportMethod coerce
#'
#' @seealso
#' - [Seurat::CreateSeuratObject()].
#' - [Seurat::Convert()].
#'
#' @examples
#' data(seurat_small)
#' data(sce_small, package = "basejump")
#'
#' ## SingleCellExperiment to seurat ====
#' object <- sce_small
#' print(object)
#' x <- as(object, "seurat")
#' class(x)
#' print(x)
#'
#' ## seurat to SingleCellExperiment ====
#' x <- as(seurat_small, "SingleCellExperiment")
#' print(x)
#'
#' ## seurat to RangedSummarizedExperiment ====
#' x <- as(seurat_small, "RangedSummarizedExperiment")
#' print(x)
#'
#' ## seurat to SummarizedExperiment ====
#' x <- as(seurat_small, "SummarizedExperiment")
#' print(x)
NULL



as.SingleCellExperiment <- function(x, ...) {
    UseMethod("as.SingleCellExperiment")
}



# Fast internal method, that ensures dimnames are sanitized.
# Sanitizes legacy seurat objects that contain Ensembl IDs in names.
as.SingleCellExperiment.seurat <- function(x, ...) {  # nolint
    validObject(x)
    names(rownames(x@raw.data)) <- NULL
    names(rownames(x@data)) <- NULL
    Seurat::Convert(from = x, to = "sce")
}



#' @rdname as
#' @name coerce,SingleCellExperiment,seurat-method
#'
#' @section `SingleCellExperiment` to `seurat`:
#' Interally [Seurat::CreateSeuratObject()] is called without applying any
#' additional filtering cutoffs, since we have already defined them during our
#' quality control analysis. Here we are passing the raw gene-level counts of
#' the filtered cells into a new `seurat` class object. Use
#' [convertGenesToSymbols()] to convert gene IDs to names (symbols).
setAs(
    from = "SingleCellExperiment",
    to = "seurat",
    function(from) {
        # Create the seurat object
        to <- CreateSeuratObject(
            raw.data = counts(from),
            project = "pointillism",
            # Already applied filtering cutoffs for cells and genes.
            min.cells = 0L,
            min.genes = 0L,
            # Using the default for UMI datasets.
            is.expr = 0L,
            meta.data = as.data.frame(colData(from))
        )

        # Check that the dimensions match exactly.
        assert_are_identical(
            x = dim(from),
            y = dim(slot(to, "raw.data"))
        )

        # Stash metadata and rowRanges into `misc` slot.
        metadata <- metadata(from)
        metadata[["sessionInfo"]] <- session_info()
        to@misc[["metadata"]] <- metadata
        to@misc[["rowRanges"]] <- rowRanges(from)

        to
    }
)



#' @rdname as
#' @name coerce,seurat,SingleCellExperiment-method
#'
#' @section `seurat` to `SingleCellExperiment`:
#' S4 coercion support for creating a `SingleCellExperiment` from a `seurat`
#' class object. The [Seurat FAQ page](https://satijalab.org/seurat/faq)
#' explains the `seurat` S4 class structure in detail. Internally, this method
#' improves the basic `Seurat::as.SingleCellExperiment()` S3 coercion method,
#' including the `object@scale.data` matrix, and will keep track of stashed
#' `rowRanges` and `metadata` if the `seurat` object was originally created
#' from a `SingleCellExperiment` (i.e. from the bcbioSingleCell package).
setAs(
    from = "seurat",
    to = "SingleCellExperiment",
    function(from) {
        to <- as.SingleCellExperiment(from)
        # Slot scaleData, if defined. Note that we're making the dimensions
        # match the other count matrices here.
        scaleData <- slot(from, "scale.data")
        if (!is.null(scaleData)) {
            scaleData <- scaleData[rownames(to), colnames(to)]
            assays(to)[["scaleData"]] <- scaleData
        }
        # Sanitize row and column data colnames using camel case.
        rowRanges(to) <- camel(rowRanges(from))
        colData(to) <- camel(colData(to))
        metadata(to) <- metadata(from)
        # Slot variable genes, if calculated.
        varGenes <- slot(from, "var.genes")
        metadata(to)[["varGenes"]] <- varGenes
        to
    }
)



#' @rdname as
#' @name coerce,seurat,RangedSummarizedExperiment-method
#'
#' @section `seurat` to `RangedSummarizedExperiment`:
#' S4 coercion support for creating a `RangedSummarizedExperiment` from a
#' `seurat` class object.
setAs(
    from = "seurat",
    to = "RangedSummarizedExperiment",
    function(from) {
        sce <- as(from, "SingleCellExperiment")
        rse <- as(sce, "RangedSummarizedExperiment")
        rse
    }
)



#' @rdname as
#' @name coerce,seurat,SummarizedExperiment-method
#' @section `seurat` to `SummarizedExperiment`:
#' S4 coercion support for creating a `SummarizedExperiment` from a `seurat`
#' class object.
setAs(
    from = "seurat",
    to = "SummarizedExperiment",
    function(from) {
        rse <- as(from, "RangedSummarizedExperiment")
        se <- as(rse, "SummarizedExperiment")
        se
    }
)
