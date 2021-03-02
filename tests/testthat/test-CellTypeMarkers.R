context("CellTypeMarkers")

test_that("CellTypeMarkers", {
    file <- system.file(
        "extdata/markers/cell-type/homo-sapiens.csv",
        package = "pointillism"
    )
    markers <- as(import(file), "DataFrame")
    gene2symbol <- Gene2Symbol(seurat)
    keep <- markers[["geneId"]] %in% gene2symbol[["geneId"]]
    expect_true(any(keep))
    markers <- markers[keep, , drop = FALSE]
    x <- CellTypeMarkers(
        object = markers,
        gene2symbol = gene2symbol
    )
    expect_s4_class(x, "CellTypeMarkers")
})
