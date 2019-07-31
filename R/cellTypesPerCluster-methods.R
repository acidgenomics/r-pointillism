#' @name cellTypesPerCluster
#' @inherit bioverbs::cellTypesPerCluster
#' @note Updated 2019-07-31.
#'
#' @inheritParams basejump::params
#' @param min `integer(1)`.
#'   Minimum number of marker genes per cluster.
#' @param max `integer(1)`.
#'   Maximum number of marker genes per cluster.
#' @param ... Additional arguments.
#'
#' @return `grouped_df`.
#' Grouped by `cluster` column, containing the count (`n`) of significant  known
#' makers per cell type.
#'
#' @examples
#' data(cellTypeMarkersList, seuratAllMarkers)
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



cellTypesPerCluster.KnownMarkers <-  # nolint
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

        ## Note that the order is important here.
        groupCols <- c("cluster", "cellType")

        data <- object %>%
            as_tibble() %>%
            ungroup() %>%
            select(!!!syms(groupCols), everything()) %>%
            mutate_at(groupCols, as.factor) %>%
            group_by(!!!syms(groupCols)) %>%
            arrange(!!sym("padj"), .by_group = TRUE) %>%
            ## Only positive markers are informative and should be used.
            filter(!!sym("avgLogFC") > 0L) %>%
            ## Use `toString()` instead of `aggregate()` for R Markdown tables.
            ## Genes are arranged by P value.
            summarize(
                n = n(),
                name = toString(!!sym("name")),
                geneID = toString(!!sym("geneID")),
                geneName = toString(!!sym("geneName"))
            ) %>%
            group_by(!!sym("cluster")) %>%
            arrange(desc(!!sym("n")), .by_group = TRUE)

        ## Apply minimum and maximum gene cutoffs.
        if (is.numeric(min) && min > 1L) {
            data <- filter(data, !!sym("n") >= !!min)
        }
        if (is.numeric(max) && max > 1L) {
            data <- filter(data, !!sym("n") <= !!max)
        }
        assert(hasRows(data))

        data
    }



#' @rdname cellTypesPerCluster
#' @export
setMethod(
    f = "cellTypesPerCluster",
    signature = signature("KnownMarkers"),
    definition = cellTypesPerCluster.KnownMarkers
)
