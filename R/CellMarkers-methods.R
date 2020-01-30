#' Cell markers
#'
#' @name CellMarkers
#' @note Updated 2020-01-30.
#'
#' @inheritParams acidroxygen::params
#'
#' @return Markers object.
NULL



## Updated 2019-08-29.
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
    x <- camelCase(x)
    x <- x[, c(group, "geneID")]
    x <- x[complete.cases(x), , drop = FALSE]
    x <- unique(x)
    ## Warn user about markers that aren't present in the gene2symbol. This is
    ## useful for informing about putative markers that aren't expressed.
    setdiff <- setdiff(x[["geneID"]], gene2symbol[["geneID"]])
    if (hasLength(setdiff)) {
        cli_alert_warning(sprintf(
            "Markers missing from gene2symbol: %s.",
            toString(setdiff, width = 200L)
        ))
    }
    intersect <- intersect(x[["geneID"]], gene2symbol[["geneID"]])
    assert(hasLength(intersect))
    keep <- x[["geneID"]] %in% intersect
    x <- x[keep, , drop = FALSE]
    x <- leftJoin(x, gene2symbol, by = "geneID")
    x <- x[, sort(colnames(x)), drop = FALSE]
    x <- x[order(x[[group]], x[["geneName"]]), , drop = FALSE]
    x <- mutateIf(x, is.character, as.factor)
    x <- split(x, f = x[[group]])
    x <- snakeCase(x)
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
    cols <- c("cellType", "geneID")
    df <- x[, cols, drop = FALSE]
    ## Generate the grouping factor necessary to perform split.
    f <- .group(df)
    split <- split(df, f = f)
    n <- vapply(X = split, FUN = nrow, FUN.VALUE = integer(1L))
    which <- which(n >= n)
    genes <- split[, "geneID"][which]
    genes <- unlist(genes, use.names = FALSE)
    genes <- sort(unique(genes))
    if (hasLength(genes)) {
        cli_alert_warning(sprintf(
            "Filtering promiscuous marker genes: %s.",
            toString(genes, width = 100L)
        ))
        keep <- !(x[["geneID"]] %in% genes)
        x <- x[keep, , drop = FALSE]
    }
    x
}



#' @describeIn CellMarkers Cell-cycle markers.
#' @export
CellCycleMarkers <-  # nolint
    function(object, gene2symbol) {
        class <- "CellCycleMarkers"
        data <- .CellMarkers(
            object = object,
            gene2symbol = gene2symbol,
            class = class
        )
        new(Class = class, data)
    }



#' @describeIn CellMarkers Cell-type markers.
#' @export
CellTypeMarkers <-  # nolint
    function(object, gene2symbol) {
        class <- "CellTypeMarkers"
        data <- .CellMarkers(
            object = object,
            gene2symbol = gene2symbol,
            class = class
        )
        new(Class = class, data)
    }
