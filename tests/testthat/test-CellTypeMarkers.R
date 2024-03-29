test_that("CellTypeMarkers", {
    gene2symbol <- GeneToSymbol(objs[["Seurat"]])
    file <- system.file(
        "extdata/markers/cell-type/homo-sapiens.csv",
        package = "AcidSingleCell"
    )
    markers <- as(import(file), "DataFrame")
    colnames(markers) <- camelCase(colnames(markers), strict = TRUE)
    keep <- markers[["geneId"]] %in% gene2symbol[["geneId"]]
    expect_true(any(keep))
    markers <- markers[keep, , drop = FALSE]
    x <- CellTypeMarkers(
        object = markers,
        gene2symbol = gene2symbol
    )
    expect_s4_class(x, "CellTypeMarkers")
})
