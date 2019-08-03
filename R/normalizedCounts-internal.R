#' Normalized counts
#'
#' Size factor and log2 normalized counts.
#'
#' @seealso
#' - `SingleCellExperiment::logcounts()`.
#' - `monocle3::normalized_counts()`.
#' - https://satijalab.org/seurat/essential_commands.html
#'
#' @noRd
.normalizedCounts <- function(object) {
    if (is(object, "Seurat")) {
        x <- object[["RNA"]]@data
    } else if (is(object, "cell_data_set")) {
        x <- monocle3::normalized_counts(
            cds = object,
            norm_method = c("log"),
            pseudocount = 1L
        )
    } else if (is(object, "SingleCellExperiment")) {
        x <- logcounts(object)
    } else {
        stop(sprintf(
            fmt = "%s is not supported.",
            class(object)[[1L]]
        ))
    }
    assert(
        is(x, "sparseMatrix"),
        hasLength(x)
    )
    x
}
