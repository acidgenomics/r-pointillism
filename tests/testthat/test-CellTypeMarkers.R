test_that("CellTypeMarkers", {
    geneToSymbol <- GeneToSymbol(objs[["Seurat"]])
    file <- system.file(
        "extdata/markers/cell-type/homo-sapiens.csv",
        package = "AcidSingleCell"
    )
    markers <- as(import(file), "DataFrame")
    colnames(markers) <- camelCase(colnames(markers), strict = TRUE)
    keep <- markers[["geneId"]] %in% geneToSymbol[["geneId"]]
    expect_true(any(keep))
    markers <- markers[keep, , drop = FALSE]
    x <- CellTypeMarkers(
        object = markers,
        geneToSymbol = geneToSymbol
    )
    expect_s4_class(x, "CellTypeMarkers")
})
