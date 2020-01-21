#' @name topMarkers
#' @inherit acidgenerics::topMarkers
#' @note Updated 2019-09-03.
#'
#' @inheritParams acidroxygen::params
#' @param n `integer(1)`.
#'   Number of genes per cluster.
#' @param direction `character(1)`.
#'   Whether to include upregulated (`"up"`; positive LFC), downregulated
#'   (`"down"`; negative LFC) or `"both"` directions of association per cluster.
#' @param ... Additional arguments.
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
#' @importFrom acidgenerics topMarkers
#' @usage topMarkers(object, ...)
#' @export
NULL



## Note that the validity method checks for sorting by adjusted P value.
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
        msg <- switch(
            EXPR = direction,
            "both" = "Including both up- and down-regulated markers.",
            "down" = "Including downregulated markers.",
            "up" = "Including upregulated markers."
        )
        message(msg)
        x <- object
        x <- SplitDataFrameList(lapply(
            X = x,
            FUN = function(x) {
                ## Subset to positive or negative correlation, if desired.
                if (identical(direction, "up")) {
                    keep <- x[["avgLogFC"]] > 0L
                    x <- x[keep, , drop = FALSE]
                } else if (identical(direction, "down")) {
                    keep <- x[["avgLogFC"]] < 0L
                    x <- x[keep, , drop = FALSE]
                }
                x <- head(x, n = n)
                x
            }
        ))
        x <- unlist(x, recursive = FALSE, use.names = FALSE)
        ranges <- x[["ranges"]]
        x[["ranges"]] <- NULL
        x[["geneID"]] <- as.character(mcols(ranges)[["geneID"]])
        x[["geneName"]] <- as.character(mcols(ranges)[["geneName"]])
        x <- mutateIf(x, is.character, as.factor)
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
