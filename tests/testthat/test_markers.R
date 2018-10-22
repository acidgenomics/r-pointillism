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



# SeuratMarkers ================================================================
test_that("SeuratMarkers", {
    ranges <- rowRanges(seurat_small)
    invisible(capture.output(
        markers <- Seurat::FindMarkers(
            object = seurat_small,
            ident.1 = "1",
            ident.2 = NULL
        )
    ))
    x <- SeuratMarkers(object = markers, ranges = ranges)
    expect_s4_class(x, "SeuratMarkers")
})

test_that("SeuratMarkersPerCluster", {
    ranges <- rowRanges(seurat_small)
    invisible(capture.output(
        markers <- Seurat::FindAllMarkers(seurat_small)
    ))
    x <- SeuratMarkersPerCluster(object = markers, ranges = ranges)
    expect_s4_class(x, "SeuratMarkersPerCluster")
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
