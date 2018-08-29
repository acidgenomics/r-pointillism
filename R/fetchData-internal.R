.fetchGeneData <- function(
    object,
    genes,
    assay = "logcounts",
    gene2symbol = FALSE,
    interestingGroups = "ident"
) {
    # Allow factor input, but coerce.
    if (is.factor(genes)) {
        genes <- as.character(genes)
    }
    assert_is_character(genes)
    assert_has_no_duplicates(genes)
    assert_is_a_string(assay)
    assert_is_subset(assay, assayNames(object))
    assert_is_a_bool(gene2symbol)
    assertFormalInterestingGroups(object, interestingGroups)

    counts <- assays(object)[[assay]]
    assert_is_subset(genes, rownames(object))
    counts <- counts[genes, , drop = FALSE]
    counts <- as.matrix(counts)

    # Convert gene IDs to gene names (symbols)
    if (isTRUE(gene2symbol) && isTRUE(.useGene2symbol(object))) {
        g2s <- gene2symbol(object)
        assertIsGene2symbol(g2s)
        g2s <- g2s[genes, , drop = FALSE]
        genes <- make.unique(g2s[["geneName"]])
        assert_are_identical(rownames(counts), g2s[["geneID"]])
        rownames(counts) <- make.unique(g2s[["geneName"]])
    }

    data <- t(counts)

    if (is.character(interestingGroups)) {
        # Always include "ident" and "sampleName" at this step.
        intgroup <- unique(c("ident", "sampleName", interestingGroups))
        intgroupData <- colData(object) %>%
            .[, intgroup, drop = FALSE] %>%
            as.data.frame()
        assert_are_identical(
            x = rownames(data),
            y = rownames(intgroupData)
        )
        data <- data %>%
            as.data.frame() %>%
            cbind(intgroupData) %>%
            as_tibble() %>%
            rownames_to_column() %>%
            uniteInterestingGroups(interestingGroups) %>%
            gather(
                key = "gene",
                value = !!sym(assay),
                !!genes
            ) %>%
            group_by(!!sym("gene"))
    }

    data
}



.fetchReducedDimData <- function(
    object,
    reducedDim,
    dimsUse = c(1L, 2L)
) {
    object <- as(object, "SingleCellExperiment")
    .assertHasIdent(object)
    assert_is_a_string(reducedDim)
    assertIsImplicitInteger(dimsUse)
    assert_is_of_length(dimsUse, 2L)

    # Reduced dimension coordinates.
    reducedDimData <- slot(object, "reducedDims")[[reducedDim]]
    if (!is.matrix(reducedDimData)) {
        stop(
            paste(reducedDim, "reduced dimension not calculated"),
            call. = FALSE
        )
    }
    reducedDimData <- camel(as.data.frame(reducedDimData))

    # Cellular barcode metrics.
    metrics <- camel(metrics(object))

    # Assert checks to make sure the cbind operation works.
    assert_are_identical(
        x = rownames(reducedDimData),
        y = rownames(metrics)
    )
    assert_are_disjoint_sets(
        x = colnames(reducedDimData),
        y = colnames(metrics)
    )

    dimCols <- colnames(reducedDimData)[dimsUse]
    assert_is_character(dimCols)

    cbind(reducedDimData, metrics) %>%
        rownames_to_column() %>%
        # Group by ident here for center calculations
        group_by(!!sym("ident")) %>%
        mutate(
            x = !!sym(dimCols[[1L]]),
            y = !!sym(dimCols[[2L]]),
            centerX = median(!!sym(dimCols[[1L]])),
            centerY = median(!!sym(dimCols[[2L]]))
        ) %>%
        ungroup() %>%
        as.data.frame() %>%
        # Ensure all columns are camel case, for consistency.
        camel() %>%
        column_to_rownames()
}



.fetchReducedDimExpressionData <- function(
    object,
    genes,
    reducedDim
) {
    assert_is_subset(genes, rownames(object))

    # Log counts
    geneData <- .fetchGeneData(
        object = object,
        genes = genes,
        assay = "logcounts"
    )

    # Expression columns
    mean <- rowMeans(geneData)
    median <- rowMedians(geneData)
    sum <- rowSums(geneData)

    # Reduced dim data
    reducedDimData <- .fetchReducedDimData(
        object = object,
        reducedDim = reducedDim
    )

    cbind(reducedDimData, mean, median, sum)
}
