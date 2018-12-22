#' @name topMarkers
#' @include globals.R
#' @inherit bioverbs::topMarkers
#' @inheritParams basejump::params
#'
#' @param n `integer(1)`.
#'   Number of genes per cluster.
#' @param direction `character(1)`.
#'   Whether to include upregulated (`"up"`; positive LFC), downregulated
#'   (`"down"`; negative LFC) or `"both"` directions of association per cluster.
#'
#' @seealso
#' - `dplyr::slice`.
#' - `dplyr::top_n`.
#'
#' @examples
#' data(all_markers_small)
#' x <- topMarkers(all_markers_small, n = 2L)
#' print(x)
NULL



#' @importFrom bioverbs topMarkers
#' @aliases NULL
#' @export
bioverbs::topMarkers



topMarkers.SeuratMarkersPerCluster <-  # nolint
    function(
        object,
        n = 10L,
        direction,
        return = c("tbl_df", "DataFrame", "SplitDataFrameList")
    ) {
        validObject(object)
        assert(isInt(n))
        direction <- match.arg(direction)
        return <- match.arg(return)

        # Using `Markers` to `tbl_df` coercion method.
        data <- as(object, "tbl_df")

        # Subset to positive or negative correlation, if desired ("direction")
        # Note that `avgDiff` has been renamed to `avgLogFC` in Seurat v2.1.
        if (direction == "up") {
            message("Including upregulated markers.")
            data <- filter(data, !!sym("avgLogFC") > 0L)
        } else if (direction == "down") {
            message("Including upregulated markers.")
            data <- filter(data, !!sym("avgLogFC") < 0L)
        } else {
            message("Including both up- and down-regulated markers.")
        }

        data <- data %>%
            # Arrange by adjusted P value.
            arrange(!!sym("padj"), .by_group = TRUE) %>%
            # Take the top rows by using slice.
            dplyr::slice(1L:n)

        # Return.
        if (return == "tbl_df") {
            message("Returning as tibble.")
            data
        } else if (return == "DataFrame") {
            message("Returning as DataFrame.")
            as(data, "DataFrame")
        } else if (return == "SplitDataFrameList") {
            message("Returning as DataFrameList.")
            data <- as(data, "DataFrame")
            assert(isSubset("cluster", colnames(data)))
            data <- split(x = data, f = data[["cluster"]])
            names(data) <- paste0("cluster", names(data))
            assert(is(data, "SplitDataFrameList"))
            data
        }
    }

formals(topMarkers.SeuratMarkersPerCluster)[["direction"]] <- direction



#' @rdname topMarkers
#' @export
setMethod(
    f = "topMarkers",
    signature = signature("SeuratMarkersPerCluster"),
    definition = topMarkers.SeuratMarkersPerCluster
)
