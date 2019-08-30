#' Generator functions
#' @include AllGenerics.R
#' @noRd
NULL



## Updated 2019-08-29.
.cellMarkers <- function(
    object,
    gene2symbol,
    class = c("CellCycleMarkers", "CellTypeMarkers")
) {
    assert(
        is(object, "DataFrame"),
        is(gene2symbol, "Gene2Symbol")
    )
    class <- match.arg(class)
    group <- switch(
        EXPR = class,
        "CellCycleMarkers" = "phase",
        "CellTypeMarkers" = "cellType"
    )
    x <- unlist(object, recursive = FALSE, use.names = FALSE)
    x <- camelCase(x)
    x <- x[, c(group, "geneID")]
    x <- x[complete.cases(x), , drop = FALSE]
    x <- unique(x)
    ## Warn user about markers that aren't present in the gene2symbol. This is
    ## useful for informing about putative markers that aren't expressed.
    setdiff <- setdiff(x[["geneID"]], gene2symbol[["geneID"]])
    if (hasLength(setdiff)) {
        stop(sprintf(
            "Markers missing from gene2symbol: %s.",
            toString(setdiff, width = 200L)
        ))
    }
    intersect <- intersect(x[["geneID"]], gene2symbol[["geneID"]])
    assert(isNonEmpty(intersect))
    keep <- x[["geneID"]] %in% intersect
    x <- x[keep, , drop = FALSE]
    x <- leftJoin(x, gene2symbol, by = "geneID")
    x <- x[order(x[[group]], x[["geneName"]]), , drop = FALSE]
    x <- split(x, f = x[[group]])
    x <- snakeCase(x)
    metadata(x) <- metadata(gene2symbol)
    x
}



#' Cell-cycle markers
#' @inheritParams acidroxygen::params
#' @note Updated 2019-07-31.
#' @export
CellCycleMarkers <-  # nolint
    function(object, gene2symbol) {
        class <- "CellCycleMarkers"
        data <- .cellMarkers(
            object = object,
            gene2symbol = gene2symbol,
            class = class
        )
        new(Class = class, data)
    }



#' Cell-type markers
#' @inheritParams acidroxygen::params
#' @note Updated 2019-07-31.
#' @export
CellTypeMarkers <-  # nolint
    function(object, gene2symbol) {
        class <- "CellTypeMarkers"
        data <- .cellMarkers(
            object = object,
            gene2symbol = gene2symbol,
            class = class
        )
        new(Class = class, data)
    }



#' Known markers
#'
#' @name KnownMarkers
#' @note Both the `markers` and `known` objects must contain Ensembl gene
#'   identifiers in the `geneID` column. We must avoid any matching operations
#'   based on the gene names, since these change often and can mismatch
#'   easily.
#' @note Updated 2019-08-07.
#'
#' @inheritParams acidroxygen::params
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
#' data(cellTypeMarkersList, seuratAllMarkers)
#'
#' ## SeuratMarkersPerCluster ====
#' markers <- seuratAllMarkers
#' known <- cellTypeMarkersList[["homoSapiens"]]
#'
#' x <- KnownMarkers(
#'     markers = markers,
#'     known = known
#' )
#' summary(x)
NULL



## FIXME Remove dplyr code.
## Updated 2019-07-31.
`KnownMarkers,SeuratMarkersPerCluster` <-  # nolint
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

        ## Coerce data.
        markers <- as(markers, "tbl_df")
        known <- as(known, "tbl_df")
        known[["geneName"]] <- NULL

        ## Determine where the known markers are located in the markers data.
        ## Here we have slotted the gene IDs inside a "ranges" column.
        assert(areIntersectingSets(markers[["geneID"]], known[["geneID"]]))
        keep <- markers[["geneID"]] %in% known[["geneID"]]
        data <- markers[keep, , drop = FALSE]

        ## Apply our alpha level cutoff.
        data <- filter(data, !!sym("padj") < !!alpha)

        ## Add the `cellType` column.
        data <- leftJoin(x = data, y = known, by = "geneID")

        ## Filter out promiscuous markers present in multiple clusters.
        ## FIXME Remove pipe here.
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
                message(sprintf(
                    "Removing promiscuous markers: %s.",
                    toString(promiscuous, width = 100L)
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
    definition = `KnownMarkers,SeuratMarkersPerCluster`
)



#' Sanitize Seurat markers
#'
#' This generator function is designed to take the original return from a Seurat
#' marker analysis and add corresponding gene annotations.
#'
#' @name SeuratMarkers
#' @note For [Seurat::FindAllMarkers()] return, rownames are correctly returned
#'   in the `gene` column.
#' @note Updated 2019-08-06.
#'
#' @inheritParams acidroxygen::params
#'
#' @param object
#'   Unmodified Seurat marker return `data.frame`.
#'   - `SeuratMarkers()`: [Seurat::FindMarkers()].
#'   - `SeuratMarkersPerCluster()`: [Seurat::FindAllMarkers()].
#' @param ranges `GRanges`.
#'   Gene annotations. Names must correspond to the rownames. The function will
#'   automatically subset the ranges and arrange them alphabetically.
#'
#' @examples
#' data(Seurat, package = "acidtest")
#'
#' ## Seurat ====
#' object <- Seurat
#' ranges <- rowRanges(object)
#'
#' ## `FindMarkers()` return.
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
#' ## `FindAllMarkers()` return.
#' invisible(capture.output(suppressWarnings(
#'     markers <- Seurat::FindAllMarkers(object)
#' )))
#' x <- SeuratMarkersPerCluster(object = markers, ranges = ranges)
#' summary(x)
NULL



## FIXME Remove dplyr code.
## Updated 2019-08-06.
`SeuratMarkers,data.frame` <-  # nolint
    function(
        object,
        ranges,
        alpha = 0.05
    ) {
        assert(
            hasRows(object),
            hasRownames(object),
            is(ranges, "GRanges"),
            isSubset(
                x = c("geneID", "geneName"),
                y = colnames(mcols(ranges))
            ),
            isAlpha(alpha)
        )

        ## Detect function from column names -----------------------------------
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
        message(sprintf("'%s()' return detected.", fun))

        ## Sanitize markers ----------------------------------------------------
        ## Coerce to tibble.
        data <- as_tibble(object, rownames = "rowname")
        ## Standardize with camel case.
        data <- camelCase(data)

        ## Seurat mode ---------------------------------------------------------
        ## Map the Seurat matrix rownames to `rownames` column in tibble.
        if (grep("^Seurat", fun)) {
            if (fun == "Seurat::FindMarkers") {
                data <- rename(data, name = !!sym("rowname"))
            } else if (fun == "Seurat::FindAllMarkers") {
                data <- data %>%
                    mutate(rowname = NULL) %>%
                    rename(name = !!sym("gene"))
            }

            ## Update legacy columns.
            if ("avgDiff" %in% colnames(data)) {
                message(
                    "Renaming legacy 'avgDiff' column to 'avgLogFC' ",
                    "(changed in Seurat v2.1)."
                )
                data[["avgLogFC"]] <- data[["avgDiff"]]
                data[["avgDiff"]] <- NULL
            }

            ## Rename P value columns to match DESeq2 conventions.
            if ("pVal" %in% colnames(data)) {
                data[["pvalue"]] <- data[["pVal"]]
                data[["pVal"]] <- NULL
            }
            if ("pValAdj" %in% colnames(data)) {
                data[["padj"]] <- data[["pValAdj"]]
                data[["pValAdj"]] <- NULL
            }
        }

        ## Ensure that required columns are present.
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
            ## `cluster` is only present in `FindAllMarkers() return`.
            data <- data %>%
                select(!!!syms(c("cluster", "name")), everything()) %>%
                group_by(!!sym("cluster")) %>%
                arrange(!!sym("padj"), .by_group = TRUE)
        } else {
            data <- data %>%
                select(!!sym("name"), everything()) %>%
                arrange(!!sym("padj"))
        }

        ## Bind ranges as column -----------------------------------------------
        data <- as(data, "DataFrame")
        ## Require that all of the markers are defined in ranges.
        assert(isSubset(unique(data[["name"]]), names(ranges)))
        data[["ranges"]] <- ranges[data[["name"]]]

        ## Add metadata and return ---------------------------------------------
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
setMethod(
    f = "SeuratMarkers",
    signature = signature("data.frame"),
    definition = `SeuratMarkers,data.frame`
)



## Updated 2019-08-06.
`SeuratMarkersPerCluster,data.frame` <-  # nolint
    `SeuratMarkers,data.frame`



#' @rdname SeuratMarkers
#' @export
setMethod(
    f = "SeuratMarkersPerCluster",
    signature = signature("data.frame"),
    definition = `SeuratMarkersPerCluster,data.frame`
)
