#' @inherit KnownMarkers-class
#' @name KnownMarkers
#'
#' @note Both the `markers` and `known` objects must contain Ensembl gene
#'   identifiers in the `geneID` column. We must avoid any matching operations
#'   based on the gene names, since these change often and can mismatch
#'   easily.
#'
#' @inheritParams basejump::params
#' @param markers `SeuratMarkers` or `SeuratMarkersPerCluster`.
#' @param known `CellTypeMarkers`. Grouped by `cellType` column. Known markers
#'   `data.frame` imported by `readCellTypeMarkers` or pulled from internal
#'   cell cycle markers data.
#' @param promiscuousThreshold `integer(1)`. Minimum number of clusters
#'   required to consider a gene marker promiscuous. Set to `0` to disable
#'   promiscuous marker filtering.
#'
#' @return `KnownMarkers`.
#'
#' @examples
#' data(all_markers_small, cell_type_markers)
#' x <- KnownMarkers(
#'     markers = all_markers_small,
#'     known = cell_type_markers$homoSapiens
#' )
#' summary(x)
NULL



KnownMarkers.SeuratMarkersPerCluster <-  # nolint
    function(
        markers,
        known,
        promiscuousThreshold = 0L
    ) {
        validObject(markers)
        validObject(known)
        assert(
            isInt(promiscuousThreshold),
            allAreNonNegative(promiscuousThreshold)
        )
        promiscuousThreshold <- as.integer(promiscuousThreshold)

        alpha <- metadata(markers)[["alpha"]]
        assert(isAlpha(alpha))

        # Coerce data.
        markers <- as(markers, "tbl_df")
        known <- as(known, "tbl_df")
        known[["geneName"]] <- NULL

        # Determine where the known markers are located in the markers data.
        # Here we have slotted the gene IDs inside a "ranges" column.
        assert(areIntersectingSets(markers[["geneID"]], known[["geneID"]]))
        keep <- markers[["geneID"]] %in% known[["geneID"]]
        data <- markers[keep, , drop = FALSE]

        # Apply our alpha level cutoff.
        data <- filter(data, !!sym("padj") < !!alpha)

        # Add the `cellType` column.
        data <- left_join(x = data, y = known, by = "geneID")

        # Filter out promiscuous markers present in multiple clusters.
        if (promiscuousThreshold > 1L) {
            cols <- c("cellType", "geneID")
            promiscuous <- data[, cols] %>%
                as_tibble() %>%
                ungroup() %>%
                group_by(!!!syms(cols)) %>%
                summarize(n = n()) %>%
                filter(!!sym("n") >= !!promiscuousThreshold) %>%
                pull("geneID")
            if (length(promiscuous)) {
                message(paste(
                    "Removing promiscuous markers:", toString(promiscuous)
                ))
                keep <- !data[["geneID"]] %in% promiscuous
                data <- data[keep, , drop = FALSE]
            }
        }

        data <- as(data, "DataFrame")
        metadata(data) <- list(
            alpha = alpha,
            version = packageVersion("pointillism"),
            date = Sys.Date()
        )
        new(Class = "KnownMarkers", data)
    }



#' @rdname KnownMarkers
#' @export
setMethod(
    f = "KnownMarkers",
    signature = signature(
        markers = "SeuratMarkersPerCluster",
        known = "CellTypeMarkers"
    ),
    definition = KnownMarkers.SeuratMarkersPerCluster
)
