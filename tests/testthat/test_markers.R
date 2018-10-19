context("Marker Analysis")

data(seurat_small, envir = environment())



# CellTypeMarkers ==============================================================
test_that("CellTypeMarkers : Mus musculus", {
    file <- system.file(
        "extdata/cell_type_markers.csv",
        package = "pointillism"
    )
    g2s <- Gene2Symbol(seurat_small)
    x <- readCellTypeMarkers(
        file = file,
        gene2symbol = g2s
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
