#' Methods for Coercing an Object to a Class
#'
#' @name coerce
#' @aliases as
#' @family S4 Object
#' @author Michael Steinbaugh
#' @importFrom methods coerce
#' @exportMethod coerce
#'
#' @return Object of new class.
#'
#' @seealso
#' - [methods::as()].
#' - [methods::canCoerce()].
#' - [Seurat::CreateSeuratObject()].
#' - [Seurat::Convert()].
NULL



# Fast internal method, that ensures dimnames are sanitized.
.as.SingleCellExperiment.seurat <- function(object) {  # nolint
    validObject(object)
    # Sanitize legacy seurat objects that contain Ensembl IDs in names.
    names(rownames(object@raw.data)) <- NULL
    names(rownames(object@data)) <- NULL
    as.SingleCellExperiment(object)
}



#' @rdname coerce
#' @name coerce,SingleCellExperiment,seurat-method
#'
#' @section `SingleCellExperiment` to `seurat`:
#' Interally [Seurat::CreateSeuratObject()] is called without applying any
#' additional filtering cutoffs, since we have already defined them during our
#' quality control analysis. Here we are passing the raw gene-level counts of
#' the filtered cells into a new `seurat` class object. Use
#' [convertGenesToSymbols()] to convert gene IDs to names (symbols).
#'
#' @examples
#' # SingleCellExperiment to seurat ====
#' x <- as(sce_small, "seurat")
#' class(x)
#' print(x)
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
        to@misc[["metadata"]] <- metadata(from)
        to@misc[["rowRanges"]] <- rowRanges(from)

        to
    }
)



#' @rdname coerce
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
#'
#' @examples
#' # seurat to SingleCellExperiment ====
#' x <- as(seurat_small, "SingleCellExperiment")
#' print(x)
setAs(
    from = "seurat",
    to = "SingleCellExperiment",
    function(from) {
        to <- .as.SingleCellExperiment.seurat(from)
        # Slot scaleData, if defined. Note that we're making the dimensions
        # match the other count matrices here.
        scaleData <- slot(from, "scale.data")
        if (!is.null(scaleData)) {
            scaleData <- scaleData[rownames(to), colnames(to)]
            assays(to)[["scaleData"]] <- scaleData
        }
        rowRanges(to) <- rowRanges(from)
        metadata(to) <- metadata(from)
        to
    }
)
