# CellCycleMarkers =============================================================
setValidity(
    Class = "CellCycleMarkers",
    method = function(object) {
        assert_is_all_of(object, "DataFrame")
        assert_has_rows(object)
        assert_are_identical(
            x = colnames(object),
            y = c("phase", "geneID", "geneName")
        )
        assert_is_factor(object[["phase"]])
        assert_is_subset(
            x = c("version", "organism", "ensemblRelease", "date"),
            y = names(metadata(object))
        )
        TRUE
    }
)



# CellTypeMarkers ==============================================================
setValidity(
    Class = "CellTypeMarkers",
    method = function(object) {
        assert_is_all_of(object, "DataFrame")
        assert_has_rows(object)
        assert_are_identical(
            x = colnames(object),
            y = c("cellType", "geneID", "geneName")
        )
        assert_is_factor(object[["cellType"]])
        assert_is_subset(
            x = c("version", "organism", "ensemblRelease", "date"),
            y = names(metadata(object))
        )
        TRUE
    }
)



# SeuratMarkers ================================================================
setValidity(
    Class = "SeuratMarkers",
    method = function(object) {
        data <- slot(object, name = "data")
        # `FindAllMarkers()`
        if ("cluster" %in% colnames(data)) {
            clusters <- data[["cluster"]]
            assert_is_factor(clusters)
            clusters <- levels(clusters)
            # Ensure that we haven't already defined `ident`.
            assert_are_identical(
                x = clusters,
                y = as.character(seq(from = 0L, to = length(clusters) - 1L))
            )
        }

        # .assertIsKnownMarkers(object)
        # requiredCols <- c(
        #     "cellType",  # bcbio
        #     "cluster",   # Seurat
        #     "geneID",    # bcbio
        #     "geneName",  # bcbio
        #     "avgLogFC",  # Seurat v2.1
        #     "padj"       # Seurat v2.1
        # )
        # assert_is_subset(requiredCols, colnames(object))
        #
        #
        #
        # .isSanitizedMarkers <- function(
        #     object,
        #     package = "Seurat"
        # ) {
        #     package <- match.arg(package)
        #
        #     # General checks -----------------------------------------------------------
        #     if (!is(object, "DataFrame")) {
        #         return(FALSE)
        #     } else if (
        #         is.null(attr(object, "vars")) ||
        #         attr(object, "vars") != "cluster"
        #     ) {
        #         return(FALSE)
        #     } else if (!"geneID" %in% colnames(object)) {
        #         return(FALSE)
        #     }
        #
        #     # Package-specific checks --------------------------------------------------
        #     if (package == "Seurat") {
        #         # Check for `Seurat::FindAllMarkers()` return.
        #         # These columns are output in an inconsistent format, so we'll sanitize
        #         # into lowerCamelCase.
        #         seuratBlacklist <- c(
        #             "avg_diff",   # Legacy, now "avg_logFC"
        #             "avg_logFC",  # Renamed in v2.1
        #             "gene",
        #             "p_val",      # We'll rename to pvalue, matching DESeq2
        #             "p_val_adj",  # New in v2.1, we'll rename to padj, matching DESeq2
        #             "pct.1",
        #             "pct.2"
        #         )
        #         if (any(seuratBlacklist %in% colnames(object))) {
        #             return(FALSE)
        #         } else {
        #             return(TRUE)
        #         }
        #     }
        # }
        #
        #
        #

        TRUE
    }
)
