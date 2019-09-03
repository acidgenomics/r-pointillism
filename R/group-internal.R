## FIXME Consider a `groupAndSplit()` type function



#' Generate a grouping factor
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
