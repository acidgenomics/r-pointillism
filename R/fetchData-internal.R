.fetchGeneData <- function(
    object,
    genes,
    assay = "logcounts",
    metadata = FALSE
) {
    validObject(object)
    rownames <- .mapGenesToRownames(object, genes)
    assert_is_subset(rownames, rownames(object))
    assert_is_a_string(assay)
    assert_is_subset(assay, assayNames(object))
    assert_is_a_bool(metadata)

    counts <- assays(object)[[assay]]
    counts <- counts[rownames, , drop = FALSE]

    # Transpose, putting the gene rownames into the columns.
    data <- Matrix::t(counts)
    # Ensure we're not accidentally coercing the matrix to a different class.
    assert_are_identical(class(counts), class(data))

    # Early return the transposed matrix, if we don't want metadata.
    # This return is used by `.fetchReducedDimExpressionData()`.
    if (!isTRUE(metadata)) {
        return(data)
    }

    # Metadata is used by the plotting functions.
    interestingGroups <- interestingGroups(object)
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
    validObject(object)
    object <- as(object, "SingleCellExperiment")
    .assertHasIdent(object)
    assert_is_a_string(reducedDim)
    assertIsImplicitInteger(dimsUse)
    assert_is_of_length(dimsUse, 2L)

    # Reduced dimension coordinates.
    reducedDimData <- slot(object, "reducedDims")[[reducedDim]]
    assert_is_matrix(reducedDimData)
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
        # Group by ident here for center calculations.
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
        column_to_rownames() %>%
        # Return with the columns sorted.
        .[, sort(colnames(.))]
}



.fetchReducedDimExpressionData <- function(
    object,
    genes,
    reducedDim
) {
    validObject(object)
    genes <- .mapGenesToRownames(object, genes)

    # Log counts
    geneMatrix <- .fetchGeneData(
        object = object,
        genes = genes,
        assay = "logcounts",
        metadata = FALSE
    )
    assert_is_matrix(geneMatrix)
    assert_are_identical(colnames(geneMatrix), genes)

    # Expression columns
    mean <- rowMeans(geneMatrix)
    median <- rowMedians(geneMatrix)
    sum <- rowSums(geneMatrix)

    # Reduced dim data
    reducedDimData <- .fetchReducedDimData(
        object = object,
        reducedDim = reducedDim
    )

    cbind(reducedDimData, mean, median, sum)
}
