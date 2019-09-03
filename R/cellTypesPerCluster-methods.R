## FIXME Update bioverbs return.



#' @name cellTypesPerCluster
#' @inherit bioverbs::cellTypesPerCluster
#' @note Updated 2019-09-03.
#'
#' @inheritParams acidroxygen::params
#' @param min `integer(1)`.
#'   Minimum number of marker genes per cluster.
#' @param max `integer(1)`.
#'   Maximum number of marker genes per cluster.
#' @param ... Additional arguments.
#'
#' @return `DataFrame`.
#'
#' @examples
#' data(cellTypeMarkersList, seuratAllMarkers)
#'
#' ## KnownMarkers ====
#' markers <- KnownMarkers(
#'     markers = seuratAllMarkers,
#'     known = cellTypeMarkersList[["homoSapiens"]]
#' )
#' x <- cellTypesPerCluster(markers)
#' print(x)
NULL



#' @rdname cellTypesPerCluster
#' @name cellTypesPerCluster
#' @importFrom bioverbs cellTypesPerCluster
#' @usage cellTypesPerCluster(object, ...)
#' @export
NULL



## Updated 2019-09-03.
`cellTypesPerCluster,KnownMarkers` <-  # nolint
    function(
        object,
        min = 1L,
        max = Inf
    ) {
        validObject(object)
        assert(
            allArePositive(c(min, max)),
            isInt(min),
            isInt(max)
        )
        x <- as(object, "DataFrame")
        ## Only positive markers are informative here.
        keep <- x[["avgLogFC"]] > 0L
        x <- x[keep, , drop = FALSE]
        x <- x[order(x[["padj"]]), , drop = FALSE]
        vars <- c("cluster", "cellType")
        f <- .group(x[, vars])
        split <- split(x = x, f = f)
        ## Summarize the data per split.
        ## Using `toString()` instead of `aggregate()` for markdown tables.
        ## Genes are arranged by adjusted P value.
        split <- SplitDataFrameList(lapply(
            X = split,
            FUN = function(x) {
                DataFrame(
                    cluster = x[["cluster"]][[1L]],
                    cellType = x[["cellType"]][[1L]],
                    name = toString(x[["name"]]),
                    geneID = toString(x[["geneID"]]),
                    geneName = toString(x[["geneName"]]),
                    n = nrow(x)
                )
            }
        ))
        x <- unlist(split, recursive = FALSE, use.names = FALSE)
        ## Apply minimum and maximum gene cutoffs.
        if (
            isTRUE(is.numeric(min)) &&
            isTRUE(min > 1L)
        ) {
            keep <- x[["n"]] >= min
            x <- x[keep, , drop = FALSE]
        }
        if (
            isTRUE(is.numeric(max)) &&
            isTRUE(max > 1L)
        ) {
            keep <- x[["n"]] <= max
            x <- x[keep, , drop = FALSE]
        }
        assert(hasRows(x))
        x <- x[order(x[["cluster"]], -x[["n"]]), , drop = FALSE]
        x
    }



#' @rdname cellTypesPerCluster
#' @export
setMethod(
    f = "cellTypesPerCluster",
    signature = signature("KnownMarkers"),
    definition = `cellTypesPerCluster,KnownMarkers`
)
