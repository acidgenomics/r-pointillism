context("Marker Analysis")

data(seurat_small, envir = environment())



# CellTypeMarkers ==============================================================
test_that("CellTypeMarkers", {
    file <- system.file(
        "extdata/cell_type_markers.csv",
        package = "pointillism"
    )
    data <- import(file)
    g2s <- Gene2Symbol(seurat_small)
    x <- CellTypeMarkers(
        object = data,
        gene2symbol = g2s
    )
    expect_is(x, "CellTypeMarkers")
})



# sanitizeSeuratMarkers ========================================================
test_that("sanitizeSeuratMarkers", {
    ranges <- rowRanges(seurat_small)

    # FindAllMarkers
    invisible(capture.output(
        markers <- Seurat::FindAllMarkers(seurat_small)
    ))
    x <- SeuratMarkers(markers = markers, ranges = ranges)
    expect_is(x, "SeuratMarkers")

    # FindMarkers
    invisible(capture.output(
        markers <- Seurat::FindMarkers(
            seurat_small,
            ident.1 = "1",
            ident.2 = NULL
        )
    ))
    x <- SeuratMarkers(markers = markers, ranges = ranges)
    expect_is(x, "SeuratMarkers")
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
            cluster = "factor",
            geneID = "character",
            geneName = "factor",
            name = "character",
            padj = "numeric",
            pct1 = "numeric",
            pct2 = "numeric",
            pvalue = "numeric"
        )
    )

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
            expect_is(x, "tbl_df")
        }
    ))
})
