# Handles either gene IDs or gene names (symbols).
.mapGenesToRownames <- function(object, genes) {
    error <- "Failed to map genes"
    genes <- as.character(genes)

    # Return if the IDs match.
    if (all(genes %in% rownames(object))) {
        return(genes)
    }

    # Get the gene-to-symbol mappings.
    # Note that gene symbols must be unique here.
    g2s <- gene2symbol(object, unique = TRUE)
    if (is.null(g2s)) {
        stop(error)
    }

    # Require that the gene2symbol rownames are identical to the object.
    # We'll be using these later for the remapping, so this step is important.
    assert_are_identical(rownames(g2s), rownames(object))

    # Attempt to detect symbols in the `genes` input vector, and map to gene IDs
    # in the rows. If we detect this, make sure all of them match.
    if (any(genes %in% g2s[["geneName"]])) {
        assert_is_subset(genes, g2s[["geneName"]])
        # Get the corresponding gene IDs and check against rownames.
        match <- match(x = genes, table = g2s[["geneName"]])
        assert_all_are_not_na(match)
        genes <- rownames(g2s[match, , drop = FALSE])
        return(genes)
    }

    # Alternatively, the user may be passing in gene IDs for an object with
    # gene symbols as the rownames. This step isn't common but we'll support it
    # and warn the user.
    if (any(genes %in% g2s[["geneID"]])) {
        assert_is_subset(genes, g2s[["geneID"]])
        # Get the corresponding gene names and check against rownames.
        match <- match(x = genes, table = g2s[["geneID"]])
        assert_all_are_not_na(match)
        genes <- rownames(g2s[match, , drop = FALSE])
        return(genes)
    }

    stop(error)
}
