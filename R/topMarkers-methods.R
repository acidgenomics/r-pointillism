## FIXME Take out the return option?
## FIXME Improve return in bioverbs.



#' @name topMarkers
#' @inherit bioverbs::topMarkers
#' @note Updated 2019-07-31.
#'
#' @inheritParams acidroxygen::params
#' @param n `integer(1)`.
#'   Number of genes per cluster.
#' @param direction `character(1)`.
#'   Whether to include upregulated (`"up"`; positive LFC), downregulated
#'   (`"down"`; negative LFC) or `"both"` directions of association per cluster.
#' @param ... Additional arguments.
#'
#' @return `DataFrame`.
#'
#' @examples
#' data(seuratAllMarkers)
#'
#' ## SeuratMarkersPerCluster ====
#' object <- seuratAllMarkers
#' x <- topMarkers(object, n = 2L)
#' print(x)
NULL



#' @rdname topMarkers
#' @name topMarkers
#' @importFrom bioverbs topMarkers
#' @usage topMarkers(object, ...)
#' @export
NULL



## Updated 2019-09-03.
`topMarkers,SeuratMarkersPerCluster` <-  # nolint
    function(
        object,
        n = 10L,
        direction
    ) {
        validObject(object)
        assert(isInt(n))
        direction <- match.arg(direction)
        x <- object
        ## Subset the split to only include the top markers.
        ## Note that the validity method checks for sorting by adjusted P value.
        x <- SplitDataFrameList(lapply(X = x, FUN = head, n = n))
        x <- unlist(x, recursive = FALSE, use.names = FALSE)
        ## Subset to positive or negative correlation, if desired.
        if (identical(direction, "up")) {
            message("Including upregulated markers.")
            keep <- x[["avgLogFC"]] > 0L
            x <- x[keep, , drop = FALSE]
        } else if (identical(direction, "down")) {
            message("Including downregulated markers.")
            keep <- x[["avgLogFC"]] < 0L
            x <- x[keep, , drop = FALSE]
        } else {
            message("Including both up- and down-regulated markers.")
        }
        ## Extract geneID and geneName columns from ranges.
        ranges <- x[["ranges"]]
        x[["ranges"]] <- NULL
        x[["geneID"]] <- mcols(ranges)[["geneID"]]
        x[["geneName"]] <- mcols(ranges)[["geneName"]]
        x
    }

formals(`topMarkers,SeuratMarkersPerCluster`)[["direction"]] <- direction



#' @rdname topMarkers
#' @export
setMethod(
    f = "topMarkers",
    signature = signature("SeuratMarkersPerCluster"),
    definition = `topMarkers,SeuratMarkersPerCluster`
)
