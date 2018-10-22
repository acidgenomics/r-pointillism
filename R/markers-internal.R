.cellMarkers <- function(
    object,
    gene2symbol,
    class = c("CellCycleMarkers", "CellTypeMarkers")
) {
    assert_is_all_of(object, "DataFrame")
    assert_is_all_of(gene2symbol, "Gene2Symbol")
    class <- match.arg(class)

    if (class == "CellCycleMarkers") {
        group <- "phase"
    } else if (class == "CellTypeMarkers") {
        group <- "cellType"
    }

    # Coerce to tibble and sanitize.
    data <- object %>%
        as_tibble() %>%
        camel() %>%
        select(!!!syms(c(group, "geneID"))) %>%
        .[complete.cases(.), , drop = FALSE]

    # Warn user about markers that aren't present in the gene2symbol.
    # This is useful for informing about putative markers that aren't expressed.
    setdiff <- setdiff(data[["geneID"]], gene2symbol[["geneID"]])
    if (length(setdiff)) {
        stop(paste(
            "Markers missing from gene2symbol:",
            printString(setdiff),
            sep = "\n"
        ))
    }
    intersect <- intersect(
        x = data[["geneID"]],
        y = gene2symbol[["geneID"]]
    )
    assert_is_non_empty(intersect)

    data <- data %>%
        filter(!!sym("geneID") %in% !!intersect) %>%
        mutate(!!sym(group) := as.factor(!!sym(group))) %>%
        unique() %>%
        left_join(
            y = as_tibble(gene2symbol),
            by = "geneID"
        ) %>%
        group_by(!!sym(group)) %>%
        arrange(!!sym("geneName"), .by_group = TRUE) %>%
        as("DataFrame")

    out <- data %>%
        split(f = .[[group]], drop = FALSE) %>%
        snake()
    metadata(out) <- metadata(gene2symbol)
    out
}




.seuratMarkers <- function(
    object,
    ranges,
    alpha = 0.05
) {
    assert_is_data.frame(object)
    assert_has_rows(object)
    assertHasRownames(object)
    assert_is_all_of(ranges, "GRanges")
    assert_is_subset(
        x = c("geneID", "geneName"),
        y = colnames(mcols(ranges))
    )
    assert_is_a_number(alpha)
    assert_all_are_in_open_range(alpha, lower = 0L, upper = 1L)

    # Detect function from column names ----------------------------------------
    seuratMarkerCols <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")
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

    # Sanitize markers ---------------------------------------------------------
    # Coerce to tibble.
    data <- as_tibble(object, rownames = "rowname")
    # Standardize with camel case.
    data <- camel(data)

    # Seurat mode --------------------------------------------------------------
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
    assert_is_subset(requiredCols, colnames(data))

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

    # Bind ranges as column ---------------------------------------------------
    data <- as(data, "DataFrame")
    # Require that all of the markers are defined in ranges.
    assert_is_subset(unique(data[["name"]]), names(ranges))
    data[["ranges"]] <- ranges[data[["name"]]]

    # Add metadata and return --------------------------------------------------
    metadata(data) <- c(
        .prototypeMetadata,
        list(
            alpha = alpha,
            sessionInfo = session_info(include_base = TRUE)
        )
    )
    data
}
