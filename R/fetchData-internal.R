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
    assert(
        isString(assay),
        isSubset(assay, assayNames(object)),
        isFlag(metadata)
    )

    rownames <- mapGenesToRownames(object, genes)
    assert(isSubset(rownames, rownames(object)))

    counts <- assays(object) %>%
        .[[assay]] %>%
        .[rownames, , drop = FALSE]

    # Transpose, putting the gene rownames into the columns.
    if (is(counts, "sparseMatrix")) {
        t <- Matrix::t
    }
    data <- t(counts)
    # Ensure we're not accidentally coercing the matrix to a different class.
    assert(identical(class(counts), class(data)))

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
    assert(
        identical(rownames(data), rownames(colData)),
        isSubset(
            x = c("ident", "interestingGroups", "sampleName"),
            y = colnames(colData)
        )
    )

    # Bind the counts and interesting groups columns.
    assert(areDisjointSets(colnames(data), colnames(colData)))
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
    assert(isNonEmpty(g2s), hasRownames(g2s))
    g2s <- as(g2s, "tbl_df")
    data <- left_join(data, g2s, by = "rowname")

    as(data, "DataFrame")
}



.fetchReducedDimData <- function(
    object,
    reducedDim,
    dimsUse = seq_len(2L)
) {
    validObject(object)
    object <- as(object, "SingleCellExperiment")
    assert(
        .hasIdent(object),
        isScalar(reducedDim),
        isIntegerish(dimsUse),
        hasLength(dimsUse, n = 2L)
    )

    # Reduced dimension coordinates.
    assert(isSubset(reducedDim, reducedDimNames(object)))
    reducedDimData <- reducedDims(object)[[reducedDim]]
    # Coerce to DataFrame, for `cbind` call below.
    reducedDimData <- as(reducedDimData, "DataFrame")

    # Cellular barcode metrics.
    colData <- metrics(object, return = "DataFrame")
    assert(
        isSubset("ident", colnames(colData)),
        identical(
            x = rownames(reducedDimData),
            y = rownames(colData)
        ),
        areDisjointSets(
            x = colnames(reducedDimData),
            y = colnames(colData)
        )
    )

    dimCols <- colnames(reducedDimData)[dimsUse]
    assert(is.character(dimCols))

    # Bind the data frames.
    data <- cbind(reducedDimData, colData)
    assert(is(data, "DataFrame"))

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
    assert(
        hasRownames(data),
        identical(rownames(data), colnames(object))
    )
    data
}

formals(.fetchReducedDimData)[c("dimsUse", "reducedDim")] <-
    list(dimsUse = dimsUse, reducedDim = reducedDim)



.fetchReducedDimExpressionData <- function(
    object,
    genes,
    reducedDim
) {
    validObject(object)
    assert(
        is.character(genes),
        isScalar(reducedDim)
    )

    rownames <- mapGenesToRownames(object, genes = genes)

    # Transposed log counts matrix, with genes in the columns.
    geneCounts <- .fetchGeneData(
        object = object,
        genes = rownames,
        assay = "logcounts",
        metadata = FALSE
    )
    assert(identical(
        x = colnames(geneCounts),
        y = as.character(rownames)
    ))

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
    assert(is(data, "DataFrame"))
    data
}
formals(.fetchReducedDimExpressionData)[["reducedDim"]] <- reducedDim
