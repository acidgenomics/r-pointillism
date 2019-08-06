#' Fetch data
#' @include globals.R
#' @noRd
NULL



## Updated 2019-08-05.
.fetchGeneData <- function(
    object,
    value = c("logcounts", "normcounts"),
    genes,
    metadata = FALSE
) {
    validObject(object)
    assert(
        isCharacter(genes),
        isFlag(metadata)
    )
    value <- match.arg(value)

    rownames <- mapGenesToRownames(object = object, genes = genes)
    assert(isSubset(rownames, rownames(object)))

    fun <- get(value, inherits = TRUE)
    assert(is.function(fun))
    counts <- fun(object)
    counts <- counts[rownames, , drop = FALSE]

    ## Transpose, putting the gene rownames into the columns.
    if (is(counts, "sparseMatrix")) {
        t <- Matrix::t
    }
    data <- t(counts)
    ## Ensure we're not accidentally coercing the matrix to a different class.
    assert(identical(class(counts), class(data)))

    ## Early return the transposed matrix, if we don't want metadata.
    ## This return is used by `.fetchReductionExpressionData()`.
    if (!isTRUE(metadata)) {
        return(data)
    }

    ## Coerce the counts matrix to a DataFrame.
    data <- as(as.matrix(data), "DataFrame")

    ## Always include "ident" and "sampleName" using `metrics` here.
    ## This ensures `sampleName` and `interestingGroups` are always defined.
    colData <- metrics(object, return = "DataFrame")
    assert(
        identical(rownames(data), rownames(colData)),
        isSubset(
            x = c("ident", "interestingGroups", "sampleName"),
            y = colnames(colData)
        )
    )

    ## Bind the counts and interesting groups columns.
    assert(areDisjointSets(colnames(data), colnames(colData)))
    data <- cbind(data, colData)

    ## Gather into long format. Here we're putting the genes into a "rowname"
    ## column. Note that this step can attempt to sanitize gene symbols (e.g.
    ## "HLA-DRA" to "HLA.DRA") here in the `as_tibble()` call, if we're not
    ## careful.
    data <- data %>%
        as_tibble(rownames = "rowname") %>%
        gather(
            key = "rowname",
            value = !!sym(value),
            !!rownames
        ) %>%
        group_by(!!sym("rowname"))

    ## Join the geneID and geneName columns by the "rowname" column.
    g2s <- Gene2Symbol(object)
    assert(isNonEmpty(g2s), hasRownames(g2s))
    g2s <- as(g2s, "tbl_df")
    data <- left_join(data, g2s, by = "rowname")
    data[["rowname"]] <- NULL
    data <- data[, sort(colnames(data))]
    as(data, "DataFrame")
}



## Updated 2019-08-02.
.fetchReductionData <- function(
    object,
    reduction,
    dimsUse = seq_len(2L)
) {
    validObject(object)
    assert(
        .hasClusters(object),
        isScalar(reduction),
        hasLength(dimsUse, n = 2L),
        all(isIntegerish(dimsUse))
    )

    ## Get reduced dimension coordinates. Map assay position to name, which
    ## we're using below to fix the naming inconsistencies in monocle3.
    if (!isString(reduction)) {
        reduction <- reducedDimNames(object)[[reduction]]
    }
    ## This step will run through on mismatch, unless we check for error above.
    reductionData <- reducedDim(object, type = reduction)
    assert(hasLength(reductionData))
    ## Handle undefined column names here, which is currently the case with
    ## monocle3 UMAP (but not PCA) output.
    if (!hasColnames(reductionData)) {
        colnames(reductionData) <-
            paste0(reduction, seq_len(ncol(reductionData)))
    }
    ## Coercing to DataFrame, for `cbind` call below.
    reductionData <- as(reductionData, "DataFrame")

    ## Cellular barcode metrics.
    colData <- metrics(object, return = "DataFrame")
    assert(
        isSubset("ident", colnames(colData)),
        identical(
            x = rownames(reductionData),
            y = rownames(colData)
        ),
        areDisjointSets(
            x = colnames(reductionData),
            y = colnames(colData)
        )
    )

    dimCols <- colnames(reductionData)[dimsUse]
    assert(is.character(dimCols))

    ## Bind the data frames.
    data <- cbind(reductionData, colData)
    assert(is(data, "DataFrame"))

    ## Coerce to long format DataFrame.
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

formals(.fetchReductionData)[c("dimsUse", "reduction")] <-
    list(dimsUse = dimsUse, reduction = reduction)



## Updated 2019-08-02.
.fetchReductionExpressionData <- function(
    object,
    genes,
    reduction
) {
    validObject(object)
    assert(
        is.character(genes),
        isScalar(reduction)
    )

    rownames <- mapGenesToRownames(object, genes = genes)

    ## Transposed log counts matrix, with genes in the columns.
    geneCounts <- .fetchGeneData(
        object = object,
        genes = rownames,
        value = "logcounts",
        metadata = FALSE
    )
    assert(identical(
        x = colnames(geneCounts),
        y = as.character(rownames)
    ))

    ## Keep the supported operations sparse.
    if (is(geneCounts, "sparseMatrix")) {
        rowMeans <- Matrix::rowMeans
        rowSums <- Matrix::rowSums
    }

    ## Calculate the expression summary columns.
    ## Note that `rowMedians` currently isn't supported for sparse data.
    mean <- rowMeans(geneCounts)
    sum <- rowSums(geneCounts)

    ## Fetch reduced dim data.
    reductionData <- .fetchReductionData(
        object = object,
        reduction = reduction
    )

    data <- cbind(reductionData, mean, sum)
    assert(is(data, "DataFrame"))
    data
}

formals(.fetchReductionExpressionData)[["reduction"]] <- reduction



## Updated 2019-08-03.
.getSeuratStash <- function(object, name) {
    assert(
        is(object, "Seurat"),
        isString(name)
    )

    misc <- slot(object, name = "misc")

    ## Early return if the `misc` slot is `NULL`.
    if (is.null(misc)) {
        return(NULL)
    }

    ## Look first directly in `object@misc` slot.
    x <- misc[[name]]
    if (!is.null(x)) {
        return(x)
    }

    ## Next, handle legacy `bcbio` stash list inside `object@misc`.
    ## As of v0.1.3, stashing directly into `object@misc`.
    if ("bcbio" %in% names(misc)) {
        x <- misc[["bcbio"]][[name]]
        if (!is.null(x)) {
            return(x)
        }
    }

    NULL
}
