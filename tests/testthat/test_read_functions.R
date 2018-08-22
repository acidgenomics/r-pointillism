context("Read Functions")



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
