# FIXME Need to update this to use new marker classes



#' Known Markers Detected
#'
#' @note Both the `all` and `known` objects must contain Ensembl gene
#'   identifiers in the `geneID` column. We must avoid any matching operations
#'   based on the gene names, since these change often and can mismatch
#'   easily.
#'
#' @name knownMarkersDetected
#' @family Marker Functions
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#' @param all `SeuratMarkers`.
#' @param known `CellTypeMarkers`. Grouped by `cellType` column. Known markers
#'   `data.frame` imported by [readCellTypeMarkers()] or pulled from internal
#'   [cell_cycle_markers] data.
#' @param alpha `scalar numeric`. Alpha cutoff (adjusted P value; false
#'   discovery rate).
#' @param filterPromiscuous `boolean`. Remove genes with poor specificity, that
#'   are present in as many clusters as defined by `promiscuousCutoff`.
#' @param promiscuousCutoff `scalar integer`. Minimum number of clusters
#'   required to consider a gene marker promiscuous.
#'
#' @return `DataFrame`.
#'
#' @examples
#' x <- knownMarkersDetected(
#'     all = all_markers_small,
#'     known = known_markers_small
#' )
#' head(x)
NULL



.knownMarkersDetected.SeuratMarkers <-  # nolint
    function(
        all,
        known,
        alpha = 0.05,
        filterPromiscuous = FALSE,
        promiscuousCutoff = 5L
    ) {
        assert_is_all_of(all, "SeuratMarkers")
        assert_is_all_of(known, "CellTypeMarkers")
        # Check for `geneID` overlap; avoid accidental organism mismatch.
        assert_are_intersecting_sets(
            x = known[["geneID"]],
            y = mcols(all@GRanges)[["geneID"]]
        )
        assert_is_a_number(alpha)
        assert_all_are_in_left_open_range(x = alpha, lower = 0L, upper = 1L)
        assert_is_a_bool(filterPromiscuous)
        # Check for valid promiscuous cutoff.
        assertIsAnImplicitInteger(promiscuousCutoff)
        assert_all_are_in_left_open_range(promiscuousCutoff, 1L)

        # Sanitize factors to character.
        all_tbl <- all %>%
            slot(name = "data") %>%
            as("tbl_df")
        # FIXME
        g2s <- all %>%
            slot(name = "GRanges") %>%
            gene2symbol()
        # FIXME
        stop("Still in development")
            all %>%
            slot("data") %>%
            as("tbl_df") %>%
            mutate_if(is.factor, as.character) %>%
            # Ensure coercion of `AsIs` column.
            mutate(!!sym("geneID") := as.character(!!sym("geneID")))
        known_tbl <- known %>%
            as("DataFrame") %>%
            as("tbl_df") %>%
            select(!!!syms(c("cellType", "geneID"))) %>%
            mutate_if(is.factor, as.character)

        # Group by cell type and arrange by P value.
        data <- all %>%
            ungroup() %>%
            # Apply alpha cutoff, before adding cell type annotations.
            filter(!!sym("padj") < !!alpha) %>%
            left_join(known, by = "geneID") %>%
            select(
                !!!syms(c("cellType", "cluster", "geneID", "geneName")),
                everything()
            ) %>%
            mutate_at(c("cellType", "cluster"), as.factor) %>%
            filter(!is.na(!!sym("cellType"))) %>%
            group_by(!!sym("cellType")) %>%
            arrange(!!sym("padj"), .by_group = TRUE)

        if (isTRUE(filterPromiscuous)) {
            # Filter out promiscuous markers present in multiple clusters.
            promiscuous <- data %>%
                ungroup() %>%
                group_by(!!!syms(c("cellType", "geneID", "geneName"))) %>%
                summarize(n = n()) %>%
                filter(!!sym("n") >= !!promiscuousCutoff) %>%
                pull("geneID")
            if (length(promiscuous)) {
                message(paste(
                    "Promiscuous markers:", toString(promiscuous)
                ))
                data <- filter(data, !!sym("geneID") %in% !!!promiscuous)
            }
        }

        data
    }
