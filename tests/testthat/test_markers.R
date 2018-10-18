context("Marker Functions")



# knownMarkersDetected =========================================================
test_that("knownMarkersDetected", {
    x <- knownMarkersDetected(
        all = all_markers_small,
        known = known_markers_small,
        filterPromiscuous = TRUE
    )
    expect_is(x, "grouped_df")
    expect_identical(dplyr::group_vars(x), "cellType")
    expect_identical(
        lapply(x, class) %>%
            .[sort(names(.))],
        list(
            avgLogFC = "numeric",
            broadClass = "character",
            cellType = "factor",
            cluster = "factor",
            description = "character",
            end = "integer",
            geneBiotype = "character",
            geneID = "character",
            geneName = "character",
            padj = "numeric",
            pct1 = "numeric",
            pct2 = "numeric",
            pvalue = "numeric",
            rowname = "character",
            seqCoordSystem = "character",
            seqnames = "character",
            start = "integer",
            strand = "character",
            width = "integer"
        )
    )
})



# readCellTypeMarkers ==========================================================
test_that("readCellTypeMarkers : Mus musculus", {
    file <- system.file(
        "extdata/cell_type_markers.csv",
        package = "pointillism"
    )
    gene2symbol <- makeGene2symbolFromEnsembl("Homo sapiens")
    x <- readCellTypeMarkers(
        file = file,
        gene2symbol = gene2symbol
    )
    group <- dplyr::group_vars(x)
    expect_identical(group, "cellType")
    y <- tibble::tibble(
        cellType = c("cell_type_1", "cell_type_2"),
        geneID = c("ENSG00000130711", "ENSG00000109099"),
        geneName = c("PRDM12", "PMP22")
    ) %>%
        dplyr::group_by(!!rlang::sym("cellType"))
    expect_equal(x, y)
})



# sanitizeSeuratMarkers ========================================================
test_that("sanitizeSeuratMarkers", {
    # Early return on sanitized data
    expect_message(
        sanitizeSeuratMarkers(
            data = all_markers_small,
            rowRanges = rowRanges(seurat_small)
        ),
        "Markers are already sanitized"
    )

    # FindAllMarkers
    invisible(capture.output(
        all <- Seurat::FindAllMarkers(seurat_small)
    ))
    x <- sanitizeSeuratMarkers(
        data = all,
        rowRanges = rowRanges(seurat_small)
    )
    expect_is(x, "grouped_df")

    # FindMarkers
    invisible(capture.output(
        ident1 <- Seurat::FindMarkers(
            seurat_small,
            ident.1 = "1",
            ident.2 = NULL
        )
    ))
    x <- sanitizeSeuratMarkers(
        data = ident1,
        rowRanges = rowRanges(seurat_small)
    )
    expect_is(x, "data.frame")
    expect_true(tibble::has_rownames(x))
})



# topMarkers ===================================================================
test_that("topMarkers : grouped_df", {
    x <- topMarkers(all_markers_small)
    expect_is(x, "grouped_df")
    expect_identical(dplyr::group_vars(x), "cluster")
    expect_identical(
        lapply(x, class) %>%
            .[sort(names(.))],
        list(
            avgLogFC = "numeric",
            broadClass = "factor",
            cluster = "factor",
            description = "factor",
            end = "integer",
            geneBiotype = "factor",
            geneID = "character",
            geneName = "factor",
            padj = "numeric",
            pct1 = "numeric",
            pct2 = "numeric",
            pvalue = "numeric",
            rowname = "character",
            seqCoordSystem = "factor",
            seqnames = "factor",
            start = "integer",
            strand = "factor",
            width = "integer"
        )
    )

    # Coding
    x <- topMarkers(all_markers_small, coding = TRUE)
    expect_is(x, "data.frame")

    # Direction
    direction <- formals(topMarkers) %>%
        .[["direction"]] %>%
        as.character() %>%
        .[-1L]
    invisible(lapply(
        X = direction,
        FUN = function(direction) {
            x <- topMarkers(
                data = all_markers_small,
                direction = direction
            )
            expect_is(x, "data.frame")
        }
    ))
})
