context("CellTypeMarkers")

test_that("CellTypeMarkers", {
    file <- system.file(
        "extdata/markers/cell-type/homo-sapiens.csv",
        package = "pointillism"
    )
    x <- CellTypeMarkers(
        object = as(import(file), "DataFrame"),
        gene2symbol = Gene2Symbol(seurat)
    )
    expect_is(x, "CellTypeMarkers")
})
