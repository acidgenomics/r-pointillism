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
    if (.isSparseMatrix(counts)) {
        t <- Matrix::t
    }
    data <- t(counts)
    # Ensure we're not accidentally coercing the matrix to a different class.
    assert_are_identical(class(counts), class(data))

    # Early return the transposed matrix, if we don't want metadata.
    # This return is used by `.fetchReducedDimExpressionData()`.
    if (!isTRUE(metadata)) {
        return(data)
    }

    # Metadata is required for the plotting functions.
    interestingGroups <- interestingGroups(object)
    assert_is_non_empty(interestingGroups)

    # Otherwise, coerce the counts matrix to a DataFrame.
    data <- as(as.matrix(data), "DataFrame")

    # Always include "ident" and "sampleName" at this step.
    intgroup <- unique(c("ident", "sampleName", interestingGroups))
    intgroupData <- colData(object)[, intgroup, drop = FALSE]
    assert_are_identical(rownames(data), rownames(intgroupData))

    # Bind the counts and interesting groups columns.
    data <- cbind(data, intgroupData)

    # Gather into long format tibble.
    # Here we're putting the genes into a "rowname" column.
    data <- data %>%
        as("tbl_df") %>%
        uniteInterestingGroups(interestingGroups) %>%
        gather(
            key = "rowname",
            value = !!sym(assay),
            !!rownames
        ) %>%
        group_by(!!sym("rowname"))

    # Join the geneID and geneName columns by the "rowname" column.
    g2s <- gene2symbol(object)
    assert_is_non_empty(g2s)
    assertHasRownames(g2s)
    g2s <- as(g2s, "tbl_df")
    data <- left_join(data, g2s, by = "rowname")

    as(data, "DataFrame")
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
    reducedDimData <- reducedDims(object)[[reducedDim]]
    assert_is_non_empty(reducedDimData)
    reducedDimData <- camel(as(reducedDimData, "DataFrame"))

    # Cellular barcode metrics.
    metrics <- camel(as(metrics(object), "DataFrame"))

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

    # Bind the reduced dim coordinates and metrics.
    data <- cbind(reducedDimData, metrics)
    assert_is_all_of(data, "DataFrame")

    # Coerce to long format tibble.
    tbl <- data %>%
        as("tbl_df") %>%
        camel() %>%
        # Group by ident here for center calculations.
        group_by(!!sym("ident")) %>%
        mutate(
            x = !!sym(dimCols[[1L]]),
            y = !!sym(dimCols[[2L]]),
            centerX = median(!!sym(dimCols[[1L]])),
            centerY = median(!!sym(dimCols[[2L]]))
        ) %>%
        # Sort the columns alphabetically.
        .[, sort(colnames(.))]

    as(tbl, "DataFrame")
}



.fetchReducedDimExpressionData <- function(
    object,
    genes,
    reducedDim
) {
    validObject(object)
    rownames <- .mapGenesToRownames(object, genes)

    # Transposed log counts matrix, with genes in the columns.
    geneCounts <- .fetchGeneData(
        object = object,
        genes = rownames,
        assay = "logcounts",
        metadata = FALSE
    )
    assert_are_identical(colnames(geneCounts), rownames)

    # Keep the supported operations sparse.
    if (.isSparseMatrix(geneCounts)) {
        rowMeans <- Matrix::rowMeans
        rowSums <- Matrix::rowSums
    }

    # Calculate the expression summary columns.
    # Note that `rowMedians()` currently isn't supported for sparse data.
    mean <- rowMeans(geneCounts)
    sum <- rowSums(geneCounts)

    # Reduced dim data
    reducedDimData <- .fetchReducedDimData(
        object = object,
        reducedDim = reducedDim
    )

    data <- cbind(reducedDimData, mean, sum)
    assert_is_all_of(data, "DataFrame")
    data
}
