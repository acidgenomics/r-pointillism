# Handle either gene IDs or gene names (symbols) with these utilities.



.mapGenesToRownames <- function(object, genes) {
    error <- "Failed to map genes"

    # Allow factor input, but coerce.
    if (is.factor(genes)) {
        genes <- as.character(genes)
    }
    assert_is_character(genes)
    assert_is_non_empty(genes)

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
        message("Remapping gene names (symbols) to identifiers")
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
        message("Remapping gene identifiers to names (symbols)")
        assert_is_subset(genes, g2s[["geneID"]])
        # Get the corresponding gene names and check against rownames.
        match <- match(x = genes, table = g2s[["geneID"]])
        assert_all_are_not_na(match)
        genes <- rownames(g2s[match, , drop = FALSE])
        return(genes)
    }

    stop(error)
}



.mapGenesToSymbols <- function(object, genes) {
    error <- "Failed to map genes"

    # Allow factor input, but coerce.
    if (is.factor(genes)) {
        genes <- as.character(genes)
    }
    assert_is_character(genes)
    assert_is_non_empty(genes)

    # Get the gene-to-symbol mappings.
    # Note that gene symbols must be unique here.
    g2s <- gene2symbol(object, unique = TRUE)
    if (is.null(g2s)) {
        stop(error)
    }

    # Require that the gene2symbol rownames are identical to the object.
    # We'll be using these later for the remapping, so this step is important.
    assert_are_identical(rownames(g2s), rownames(object))

    # User passed in gene identifiers.
    if (any(genes %in% g2s[["geneID"]])) {
        message("Remapping gene identifiers to names (symbols)")
        assert_is_subset(genes, g2s[["geneID"]])
        # Get the corresponding gene names and check against rownames.
        match <- match(x = genes, table = g2s[["geneID"]])
        assert_all_are_not_na(match)
        genes <- g2s[match, "geneName", drop = TRUE]
        assert_is_character(genes)
        return(genes)
    }

    # User passed in gene names (symbols).
    # Attempt to detect pass in of non-unique/ambiguous symbols here, and warn
    # the user when this happens.
    if (any(genes %in% g2s[["geneName"]])) {
        message("Checking to ensure symbols are non-ambiguous")
        assert_is_subset(genes, g2s[["geneName"]])
        # Check the original un-sanitized gene symbols, for setdiff.
        cleanGenes <- mcols(rowRanges(object))[["geneName"]]
        assert_is_non_empty(cleanGenes)
        cleanGenes <- unique(as.character(cleanGenes))
        setdiff <- setdiff(genes, cleanGenes)
        if (length(setdiff)) {
            warning(paste(
                "Ambiguous genes:", toString(setdiff)
            ))
        }
        return(genes)
    }

    stop(error)
}
