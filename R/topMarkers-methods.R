#' @name topMarkers
#' @inherit AcidGenerics::topMarkers
#' @note Updated 2021-03-03.
#'
#' @inheritParams AcidRoxygen::params
#' @param direction `character(1)`.
#'   Whether to include upregulated (`"up"`; positive LFC), downregulated
#'   (`"down"`; negative LFC) or `"both"` directions of association per cluster.
#' @param n `integer(1)`.
#'   Number of genes per cluster.
#' @param ... Additional arguments.
#'
#' @examples
#' data(seuratAllMarkers)
#'
#' ## SeuratMarkersPerCluster ====
#' object <- seuratAllMarkers
#' x <- topMarkers(
#'     object = object,
#'     direction = "up",
#'     n = 2L
#' )
#' print(x)
NULL



## Note that the validity method checks for sorting by adjusted P value.
## Updated 2021-03-03.
`topMarkers,SeuratMarkersPerCluster` <-  # nolint
    function(
        object,
        direction,
        n = 10L
    ) {
        validObject(object)
        lfcCol <- "avgLog2Fc"
        assert(
            isInt(n),
            isSubset(lfcCol, colnames(object[[1L]]))
        )
        direction <- match.arg(direction)
        alertInfo(switch(
            EXPR = direction,
            "both" = "Including both up- and down-regulated markers.",
            "down" = "Including downregulated markers.",
            "up" = "Including upregulated markers."
        ))
        x <- object
        x <- SplitDataFrameList(lapply(
            X = x,
            FUN = function(x) {
                ## Subset to positive or negative correlation, if desired.
                if (identical(direction, "up")) {
                    keep <- x[[lfcCol]] > 0L
                    x <- x[keep, , drop = FALSE]
                } else if (identical(direction, "down")) {
                    keep <- x[[lfcCol]] < 0L
                    x <- x[keep, , drop = FALSE]
                }
                x <- head(x, n = n)
                x
            }
        ))
        x <- unlist(x, recursive = FALSE, use.names = FALSE)
        ranges <- x[["ranges"]]
        x[["ranges"]] <- NULL
        x[["geneId"]] <- as.character(mcols(ranges)[["geneId"]])
        x[["geneName"]] <- as.character(mcols(ranges)[["geneName"]])
        x <- mutateIf(x, is.character, as.factor)
        x
    }

args <- "direction"
formals(`topMarkers,SeuratMarkersPerCluster`)[args] <-
    .formalsList[args]
rm(args)



#' @rdname topMarkers
#' @export
setMethod(
    f = "topMarkers",
    signature = signature("SeuratMarkersPerCluster"),
    definition = `topMarkers,SeuratMarkersPerCluster`
)
