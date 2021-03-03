## Updated 2021-03-02.
.CellMarkers <- function(  # nolint
    object,
    gene2symbol,
    class = c("CellCycleMarkers", "CellTypeMarkers")
) {
    assert(
        is(object, "DataFrame"),
        is(gene2symbol, "Gene2Symbol")
    )
    class <- match.arg(class)
    group <- switch(
        EXPR = class,
        "CellCycleMarkers" = "phase",
        "CellTypeMarkers" = "cellType"
    )
    x <- object
    x <- camelCase(x, strict = TRUE)
    cols <- c(group, "geneId")
    x <- x[, cols]
    x <- x[complete.cases(x), , drop = FALSE]
    x <- unique(x)
    ## Warn user about markers that aren't present in the gene2symbol. This is
    ## useful for informing about putative markers that aren't expressed.
    setdiff <- setdiff(x[["geneId"]], gene2symbol[["geneId"]])
    if (hasLength(setdiff)) {
        alertWarning(sprintf(
            "Markers missing from gene2symbol: %s.",
            toString(setdiff, width = 200L)
        ))
    }
    intersect <- intersect(x[["geneId"]], gene2symbol[["geneId"]])
    assert(hasLength(intersect))
    keep <- x[["geneId"]] %in% intersect
    x <- x[keep, , drop = FALSE]
    x <- leftJoin(
        x = x,
        y = as(gene2symbol, "DataFrame"),
        by = "geneId"
    )
    x <- x[, sort(colnames(x)), drop = FALSE]
    x <- x[order(x[[group]], x[["geneName"]]), , drop = FALSE]
    x <- mutateIf(x, is.character, as.factor)
    x <- split(x, f = x[[group]])
    names(x) <- snakeCase(names(x))
    ## Specific fix for G2/M input (cell-cycle markers).
    names(x) <- sub("g2_slash_m", "g2m", names(x))
    ## FIXME This isn't passing through to the valiidty check sufficiently...
    ## FIXME makeGene2SymbolFromEnsembl isn't defining metadata correctly.
    ## Need to update this in AcidGenomes.
    metadata(x) <- metadata(gene2symbol)
    x
}



## Using approach in `uniteInterestingGroups()` to generate the
## factor grouping column to apply split.
## Updated 2020-01-30.
.filterPromiscuousMarkers <- function(x, n) {
    assert(
        is(x, "DataFrame"),
        isInt(n),
        isTRUE(n > 1L)
    )
    cols <- c("cellType", "geneId")
    df <- x[, cols, drop = FALSE]
    ## Generate the grouping factor necessary to perform split.
    f <- .group(df)
    split <- split(df, f = f)
    n <- vapply(X = split, FUN = nrow, FUN.VALUE = integer(1L))
    which <- which(n >= n)
    genes <- split[, "geneId"][which]
    genes <- unlist(genes, use.names = FALSE)
    genes <- sort(unique(genes))
    if (hasLength(genes)) {
        alertWarning(sprintf(
            "Filtering promiscuous marker genes: %s.",
            toString(genes, width = 100L)
        ))
        keep <- !(x[["geneId"]] %in% genes)
        x <- x[keep, , drop = FALSE]
    }
    x
}
