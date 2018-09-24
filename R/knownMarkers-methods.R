#' Known Markers Detected
#'
#' @note Both the `all` and `known` objects must contain Ensembl gene
#'   identifiers in the `geneID` column. We must avoid any matching operations
#'   based on the gene names, since these change often and can mismatch
#'   easily.
#'
#' @name knownMarkers
#' @family Marker Functions
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#' @param all `SeuratMarkers`.
#' @param known `CellTypeMarkers`. Grouped by `cellType` column. Known markers
#'   `data.frame` imported by [readCellTypeMarkers()] or pulled from internal
#'   [cell_cycle_markers] data.
#' @param promiscuousThreshold `scalar integer`. Minimum number of clusters
#'   required to consider a gene marker promiscuous. Set to `0` to disable
#'   promiscuous marker filtering.
#'
#' @return `KnownSeuratMarkers`.
#'
#' @examples
#' # FIXME Need to update this.
#' \dontrun{
#' x <- knownMarkers(
#'     all = all_markers_small,
#'     known = known_markers_small
#' )
#' head(x)
#' }
NULL



.knownMarkers.SeuratMarkers <-  # nolint
    function(
        all,
        known,
        promiscuousThreshold = 5L
    ) {
        validObject(all)
        validObject(known)
        assertIsAnImplicitInteger(promiscuousThreshold)
        assert_all_are_non_negative(promiscuousThreshold)
        promiscuousThreshold <- as.integer(promiscuousThreshold)

        alpha <- metadata(all)[["alpha"]]
        assertIsAlpha(alpha)

        # Determine where the known markers are located in the all markers data.
        # Here we have slotted the gene IDs inside a "ranges" column.
        allGenes <- all %>%
            .[["ranges"]] %>%
            mcols() %>%
            .[["geneID"]]
        knownGenes <- known[["geneID"]]
        assert_are_intersecting_sets(knownGenes, allGenes)
        keep <- allGenes %in% knownGenes
        data <- all[keep, , drop = FALSE]

        # Apply our alpha level cutoff.
        keep <- data[["padj"]] < alpha
        data <- data[keep, , drop = FALSE]

        # Add the `cellType` column.
        map <- left_join(
            x = DataFrame(
                name = data[["name"]],
                geneID = mcols(data[["ranges"]])[["geneID"]]
            ),
            y = known %>%
                as("DataFrame") %>%
                select(!!!syms(c("cellType", "geneID"))),
            by = "geneID"
        )
        data[["cellType"]] <- map[["cellType"]]

        # Filter out promiscuous markers present in multiple clusters.
        if (promiscuousThreshold > 1L) {
            cols <- c("cellType", "name")
            promiscuous <- data[, cols] %>%
                as("tbl_df") %>%
                ungroup() %>%
                group_by(!!!syms(cols)) %>%
                summarize(n = n()) %>%
                filter(!!sym("n") >= !!promiscuousThreshold) %>%
                pull("name")
            if (length(promiscuous)) {
                message(paste(
                    "Removing promiscuous markers:", toString(promiscuous)
                ))
                keep <- !data[["name"]] %in% promiscuous
                data <- data[keep, , drop = FALSE]
            }
        }

        metadata(data) <- list(
            alpha = alpha,
            version = packageVersion("pointillism"),
            date = Sys.Date()
        )
        new(Class = "KnownSeuratMarkers", data)
    }



#' @rdname knownMarkers
#' @export
setMethod(
    f = "knownMarkers",
    signature = signature(
        all = "SeuratMarkers",
        known = "CellTypeMarkers"
    ),
    definition = .knownMarkers.SeuratMarkers
)
