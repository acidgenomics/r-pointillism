#' Fetch Data
#' @include globals.R
#' @noRd
NULL



.fetchGeneData <- function(
    object,
    genes,
    assay = "logcounts",
    metadata = FALSE
) {
    validObject(object)
    object <- as(object, "SingleCellExperiment")
    assert_is_a_string(assay)
    assert_is_subset(assay, assayNames(object))
    assert_is_a_bool(metadata)

    rownames <- mapGenesToRownames(object, genes)
    assert_is_subset(rownames, rownames(object))

    counts <- assays(object) %>%
        .[[assay]] %>%
        .[rownames, , drop = FALSE]

    # Transpose, putting the gene rownames into the columns.
    if (is(counts, "sparseMatrix")) {
        t <- Matrix::t
    }
    data <- t(counts)
    # Ensure we're not accidentally coercing the matrix to a different class.
    assert_are_identical(class(counts), class(data))

    # Early return the transposed matrix, if we don't want metadata.
    # This return is used by `.fetchReducedDimExpressionData`.
    if (!isTRUE(metadata)) {
        return(data)
    }

    # Coerce the counts matrix to a DataFrame.
    data <- as(as.matrix(data), "DataFrame")

    # Always include "ident" and "sampleName" using `metrics` here.
    # This ensures `sampleName` and `interestingGroups` are always defined.
    colData <- metrics(object, return = "DataFrame")
    assert_are_identical(rownames(data), rownames(colData))
    assert_is_subset(
        x = c("ident", "interestingGroups", "sampleName"),
        y = colnames(colData)
    )

    # Bind the counts and interesting groups columns.
    assert_are_disjoint_sets(colnames(data), colnames(colData))
    data <- cbind(data, colData)

    # Gather into long format.
    # Here we're putting the genes into a "rowname" column.
    data <- data %>%
        as_tibble(rownames = "rowname") %>%
        gather(
            key = "rowname",
            value = !!sym(assay),
            !!rownames
        ) %>%
        group_by(!!sym("rowname"))

    # Join the geneID and geneName columns by the "rowname" column.
    g2s <- Gene2Symbol(object)
    assert_is_non_empty(g2s)
    assertHasRownames(g2s)
    g2s <- as(g2s, "tbl_df")
    data <- left_join(data, g2s, by = "rowname")

    as(data, "DataFrame")
}



.fetchReducedDimData <- function(
    object,
    reducedDim,
    dimsUse
) {
    validObject(object)
    object <- as(object, "SingleCellExperiment")
    .assertHasIdent(object)
    assert_is_scalar(reducedDim)
    assertIsImplicitInteger(dimsUse)
    assert_is_of_length(dimsUse, 2L)

    # Reduced dimension coordinates.
    assert_is_subset(reducedDim, reducedDimNames(object))
    reducedDimData <- reducedDims(object)[[reducedDim]]
    # Coerce to DataFrame, for `cbind` call below.
    reducedDimData <- as(reducedDimData, "DataFrame")

    # Cellular barcode metrics.
    colData <- metrics(object, return = "DataFrame")
    assert_is_subset("ident", colnames(colData))

    # Assert checks to make sure the cbind operation works.
    assert_are_identical(
        x = rownames(reducedDimData),
        y = rownames(colData)
    )
    assert_are_disjoint_sets(
        x = colnames(reducedDimData),
        y = colnames(colData)
    )

    dimCols <- colnames(reducedDimData)[dimsUse]
    assert_is_character(dimCols)

    # Bind the data frames.
    data <- cbind(reducedDimData, colData)
    assert_is_all_of(data, "DataFrame")

    # Coerce to long format DataFrame.
    data <- data %>%
        as_tibble() %>%
        group_by(!!sym("ident")) %>%
        mutate(
            x = !!sym(dimCols[[1L]]),
            y = !!sym(dimCols[[2L]]),
            centerX = median(!!sym(dimCols[[1L]])),
            centerY = median(!!sym(dimCols[[2L]]))
        ) %>%
        as("DataFrame")
    assertHasRownames(data)
    assert_are_identical(rownames(data), colnames(object))
    data
}
formals(.fetchReducedDimData)[c(
    "dimsUse",
    "reducedDim"
)] <- list(
    dimsUse = dimsUse,
    reducedDim = reducedDim
)



.fetchReducedDimExpressionData <- function(
    object,
    genes,
    reducedDim
) {
    validObject(object)
    assert_is_character(genes)
    assert_is_scalar(reducedDim)

    rownames <- mapGenesToRownames(object, genes = genes)

    # Transposed log counts matrix, with genes in the columns.
    geneCounts <- .fetchGeneData(
        object = object,
        genes = rownames,
        assay = "logcounts",
        metadata = FALSE
    )
    assert_are_identical(
        x = colnames(geneCounts),
        y = as.character(rownames)
    )

    # Keep the supported operations sparse.
    if (is(geneCounts, "sparseMatrix")) {
        rowMeans <- Matrix::rowMeans
        rowSums <- Matrix::rowSums
    }

    # Calculate the expression summary columns.
    # Note that `rowMedians` currently isn't supported for sparse data.
    mean <- rowMeans(geneCounts)
    sum <- rowSums(geneCounts)

    # Fetch reduced dim data.
    reducedDimData <- .fetchReducedDimData(
        object = object,
        reducedDim = reducedDim
    )

    data <- cbind(reducedDimData, mean, sum)
    assert_is_all_of(data, "DataFrame")
    data
}
formals(.fetchReducedDimExpressionData)[["reducedDim"]] <- reducedDim
