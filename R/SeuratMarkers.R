#' @inherit SeuratMarkers-class
#' @export
#'
#' @description
#' This generator function is designed to take the original return from a Seurat
#' marker analysis and add corresponding gene annotations.
#'
#' @note For `Seurat::FindAllMarkers` return, rownames are correctly returned
#'   in the `gene` column.
#'
#' @inheritParams basejump::params
#' @param object `data.frame`. Unmodified `Seurat::FindMarkers` or
#'   `Seurat::FindAllMarkers` return.
#' @param ranges `GRanges`. Gene annotations. Names must correspond to the
#'   rownames. The function will automatically subset the ranges and arrange
#'   them alphabetically.
#'
#' @examples
#' data(seurat_small)
#' object <- seurat_small
#' ranges <- rowRanges(object)
#'
#' ## `FindMarkers` return.
#' invisible(capture.output(
#'     markers <- Seurat::FindMarkers(
#'         object = object,
#'         ident.1 = "1",
#'         ident.2 = NULL
#'     )
#' ))
#' x <- SeuratMarkers(object = markers, ranges = ranges)
#' summary(x)
#'
#' ## `FindAllMarkers` return.
#' invisible(capture.output(
#'     markers <- Seurat::FindAllMarkers(object)
#' ))
#' x <- SeuratMarkersPerCluster(object = markers, ranges = ranges)
#' summary(x)
SeuratMarkers <-  # nolint
    function(
        object,
        ranges,
        alpha = 0.05
    ) {
        assert(
            is.data.frame(object),
            hasRows(object),
            hasRownames(object),
            is(ranges, "GRanges"),
            isSubset(
                x = c("geneID", "geneName"),
                y = colnames(mcols(ranges))
            ),
            isAlpha(alpha)
        )

        # Detect function from column names ------------------------------------
        seuratMarkerCols <-
            c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")
        if (identical(
            colnames(object),
            seuratMarkerCols
        )) {
            perCluster <- FALSE
            fun <- "Seurat::FindMarkers"
        } else if (identical(
            colnames(object),
            c(seuratMarkerCols, "cluster", "gene")
        )) {
            perCluster <- TRUE
            fun <- "Seurat::FindAllMarkers"
        }
        message(paste0("`", fun, "` return detected."))

        # Sanitize markers -----------------------------------------------------
        # Coerce to tibble.
        data <- as_tibble(object, rownames = "rowname")
        # Standardize with camel case.
        data <- camel(data)

        # Seurat mode ----------------------------------------------------------
        # Map the Seurat matrix rownames to `rownames` column in tibble.
        if (grep("^Seurat", fun)) {
            if (fun == "Seurat::FindMarkers") {
                data <- rename(data, name = !!sym("rowname"))
            } else if (fun == "Seurat::FindAllMarkers") {
                data <- data %>%
                    mutate(rowname = NULL) %>%
                    rename(name = !!sym("gene"))
            }

            # Update legacy columns.
            if ("avgDiff" %in% colnames(data)) {
                message(paste(
                    "Renaming legacy `avgDiff` column to `avgLogFC`",
                    "(changed in Seurat v2.1)."
                ))
                data[["avgLogFC"]] <- data[["avgDiff"]]
                data[["avgDiff"]] <- NULL
            }

            # Rename P value columns to match DESeq2 conventions.
            if ("pVal" %in% colnames(data)) {
                data[["pvalue"]] <- data[["pVal"]]
                data[["pVal"]] <- NULL
            }
            if ("pValAdj" %in% colnames(data)) {
                data[["padj"]] <- data[["pValAdj"]]
                data[["pValAdj"]] <- NULL
            }
        }

        # Ensure that required columns are present.
        requiredCols <- c(
            "name",
            "pct1",
            "pct2",
            "avgLogFC",     # Seurat v2.1.
            "padj",
            "pvalue"        # Renamed from `p_val`.
        )
        assert(isSubset(requiredCols, colnames(data)))

        if (isTRUE(perCluster)) {
            # `cluster` is only present in `FindAllMarkers() return`.
            data <- data %>%
                select(!!!syms(c("cluster", "name")), everything()) %>%
                group_by(!!sym("cluster")) %>%
                arrange(!!sym("padj"), .by_group = TRUE)
        } else {
            data <- data %>%
                select(!!sym("name"), everything()) %>%
                arrange(!!sym("padj"))
        }

        # Bind ranges as column ------------------------------------------------
        data <- as(data, "DataFrame")
        # Require that all of the markers are defined in ranges.
        assert(isSubset(unique(data[["name"]]), names(ranges)))
        data[["ranges"]] <- ranges[data[["name"]]]

        # Add metadata and return ----------------------------------------------
        metadata(data) <- c(
            .prototypeMetadata,
            list(
                alpha = alpha,
                sessionInfo = session_info(include_base = TRUE)
            )
        )

        if (isTRUE(perCluster)) {
            out <- split(x = data, f = data[["cluster"]], drop = FALSE)
            names(out) <- paste0("cluster", names(out))
            metadata(out) <- metadata(data)
            new(Class = "SeuratMarkersPerCluster", out)
        } else {
            rownames(data) <- data[["name"]]
            data[["name"]] <- NULL
            new(Class = "SeuratMarkers", data)
        }
    }



#' @rdname SeuratMarkers
#' @export
SeuratMarkersPerCluster <- SeuratMarkers  # nolint
