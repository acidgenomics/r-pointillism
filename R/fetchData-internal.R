#' Bind all dimension reduction matrices into a single data frame.
#' @note Updated 2019-08-06.
#' @noRd
.bindReducedDims <- function(object) {
    reducedDims <- reducedDims(object)
    assert(hasNames(reducedDims))
    ## Handle undefined column names here, which is currently the case with
    ## monocle3 UMAP (but not PCA) output.
    if (!isTRUE(all(bapply(X = reducedDims, FUN = hasColnames)))) {
        list <- mapply(
            name = names(reducedDims),
            assay = reducedDims,
            FUN = function(name, assay) {
                if (!hasColnames(assay)) {
                    colnames(assay) <- paste0(name, seq_len(ncol(assay)))
                }
                assay
            }
        )
        reducedDims <- as(list, "SimpleList")
    }
    out <- do.call(what = cbind, args = reducedDims)
    out <- as(out, "DataFrame")
    assert(
        hasValidDimnames(out),
        identical(x = colnames(object), y = rownames(out))
    )
    out
}



## Updated 2019-08-12.
.fetchGeneData <- function(
    object,
    genes,
    value = c("logcounts", "normcounts"),
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
    ## Transpose, putting genes into the columns.
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
    ## Gather into long format.
    ##
    ## Is there a base R / Bioconductor way to perform this step? If so, we can
    ## remove the tidyr dependency from the package.
    ##
    ## `reshape2::melt()` is also worth a look as an alternative.
    ##
    ## Here we're putting the genes into a "rowname" column. Note that this step
    ## can attempt to sanitize gene symbols (e.g. "HLA-DRA" to "HLA.DRA") here
    ## either during an `as.data.frame()` or `as_tibble()` call. To get around
    ## this, we're coercing the S4 DataFrame using `as()` and then coercing to a
    ## tibble later in the chain.
    data <- as(data, "data.frame")
    data[["rowname"]] <- rownames(data)

    data <- gather(
            data = data,
            key = "rowname",
            value = !!sym(value),
            !!rownames,
            factor_key = TRUE
    )
    data <- as(data, "DataFrame")
    ## Double checking here for accidental base R sanitization of symbols.
    assert(isSubset(rownames, data[["rowname"]]))
    ## Join the geneID and geneName columns by the "rowname" column.
    g2s <- Gene2Symbol(object)
    assert(isNonEmpty(g2s), hasRownames(g2s))
    g2s <- as(g2s, "DataFrame")
    g2s[["rowname"]] <- rownames(g2s)
    data <- left_join(data, g2s, by = "rowname")
    data <- data[, unique(c("rowname", sort(colnames(data))))]
    data
}



## Updated 2019-08-23.
.fetchReductionData <- function(
    object,
    reduction = 1L,
    dims = seq_len(2L)
) {
    validObject(object)
    assert(
        .hasClusters(object),
        isScalar(reduction),
        hasLength(dims, n = 2L),
        all(isIntegerish(dims))
    )
    ## Get reduced dimension coordinates. Map assay position to name, which
    ## we're using below to fix the naming inconsistencies in monocle3.
    if (!isString(reduction)) {
        reduction <- reducedDimNames(object)[[reduction]]
    }
    ## This step will run through on mismatch, unless we check for error above.
    redData <- reducedDim(object, type = reduction)
    assert(hasLength(redData))
    ## Handle undefined column names here, which is currently the case with
    ## monocle3 UMAP (but not PCA) output.
    if (!hasColnames(redData)) {
        colnames(redData) <- paste0(reduction, seq_len(ncol(redData)))
    }
    ## Coercing to DataFrame, for `cbind` call below.
    redData <- as(redData, "DataFrame")
    ## Cellular barcode metrics.
    colData <- metrics(object, return = "DataFrame")
    assert(
        isSubset("ident", colnames(colData)),
        identical(
            x = rownames(redData),
            y = rownames(colData)
        ),
        areDisjointSets(
            x = colnames(redData),
            y = colnames(colData)
        )
    )
    dimCols <- colnames(redData)[dims]
    assert(is.character(dimCols))
    ## Bind the data frames.
    data <- cbind(redData, colData)
    assert(is(data, "DataFrame"))
    ## Split by cluster.
    f <- data[["ident"]]
    names(f) <- rownames(data)
    split <- split(x = data, f = f)
    split <- SplitDataFrameList(lapply(
        X = split,
        FUN = function(x) {
            x[["x"]] <- x[[dimCols[[1L]]]]
            x[["y"]] <- x[[dimCols[[2L]]]]
            x[["centerX"]] <- median(x[["x"]])
            x[["centerY"]] <- median(x[["y"]])
            x
        }
    ))
    ## Note that this is using S4 `unsplit()` method defined in basejump.
    ## Method support for `SplitDataFrameList` needs to be improved in IRanges.
    data <- unsplit(value = split, f = f)
    assert(
        hasRownames(data),
        identical(rownames(data), colnames(object))
    )
    data
}

formals(.fetchReductionData)[c("dims", "reduction")] <-
    list(dims = dims, reduction = reduction)



## Updated 2019-08-23.
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



## Updated 2019-08-23.
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
