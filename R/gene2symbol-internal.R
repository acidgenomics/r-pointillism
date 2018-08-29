# Handles either gene IDs or gene names (symbols).
.mapGenesToRownames <- function(object, genes) {
    genes <- as.character(genes)

    # Return if the IDs match.
    if (all(genes %in% rownames(object))) {
        return(genes)
    }

    # Look to see if we need to map symbols to gene IDs.
    g2s <- gene2symbol(object)
    rr <- rowRanges(object)
    any(duplicated(mcols(rr)[["geneName"]]))
    any(duplicated(g2s[["geneName"]]))
    # FIXME Make sure gene2symbol returns unique symbols


    stop("Failed to map genes")
}



# Convert gene vector to symbols, if necessary.
# Currently used by plotting functions to use gene symbols with ggplot.
.mapGenes <- function(object, genes) {
    stopifnot(is(object, "SingleCellExperiment"))
    if (isTRUE(.useGene2symbol(object))) {
        g2s <- gene2symbol(object)
        if (length(g2s)) {
            g2s <- g2s[genes, , drop = FALSE]
            genes <- make.unique(g2s[["geneName"]])
        }
    }
    genes
}



#' Use Gene-to-Symbol Mappings
#'
#' Determine whether we should use stashed gene-to-symbol mappings.
#'
#' @noRd
#'
#' @examples
#' # TRUE
#' .useGene2symbol(sce_small)
#'
#' # FALSE
#' .useGene2symbol(seurat_small)
.useGene2symbol <- function(object) {
    geneNames <- tryCatch(gene2symbol(object)[["geneName"]])
    if (!length(geneNames)) {
        return(FALSE)
    }
    !any(geneNames %in% rownames(object))
}
