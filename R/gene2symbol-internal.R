# Convert gene vector to symbols, if necessary.
# Currently used by plotting functions to use gene symbols with ggplot.
.mapGenes <- function(object, genes) {
    stopifnot(is(object, "SingleCellExperiment"))
    if (isTRUE(.useGene2symbol(object))) {
        g2s <- gene2symbol(object)
        if (length(g2s)) {
            g2s <- g2s[genes, , drop = FALSE]
            genes <- make.unique(g2s[["geneName"]])
            stopifnot(all(genes %in% data[["gene"]]))
        }
    }
    genes
}



# Determine whether we should use stashed gene-to-symbol mappings.
.useGene2symbol <- function(object) {
    geneName <- as.character(
        suppressWarnings(
            gene2symbol(object)[["geneName"]]
        )
    )
    !any(geneName %in% rownames(object))
}
