#' Filter promiscuous markers
#'
#' @note Updated 2022-05-11.
#' @noRd
#'
#' @details
#' Using approach in `uniteInterestingGroups()` to generate the factor grouping
#' column to apply split.
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



#' Generate a grouping factor
#'
#' Use this internal function to define grouping factor for `split()` call.
#'
#' Use `[` subsetting to request specific columns.
#'
#' @noRd
#' @note Updated 2019-09-01.
.group <- function(x) {
    f <- apply(
        X = as.data.frame(x),
        MARGIN = 1L,
        FUN = paste,
        collapse = "."
    )
    f <- make.names(f, unique = FALSE)
    f <- as.factor(f)
    f
}
