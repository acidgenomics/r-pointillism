context("CellTypeMarkers")

## Error in validObject(.Object) :
##   invalid class "CellTypeMarkers" object: identical(x = lapply(object[[1L]], class), y = list(cellType = "factor",  .... is not TRUE.
## Calls: CellTypeMarkers -> new -> initialize -> initialize -> validObject
## Backtrace:
##     █
##  1. └─pointillism::CellTypeMarkers(object = markers, gene2symbol = gene2symbol)
##  2.   └─methods::new(Class = class, data) /data00/draco/acidbase/packages/pointillism/R/AllGenerators.R:114:8
##  3.     ├─methods::initialize(value, ...)
##  4.     └─methods::initialize(value, ...)
##  5.       └─methods::validObject(.Object)

test_that("CellTypeMarkers", {
    file <- system.file(
        "extdata/markers/cell-type/homo-sapiens.csv",
        package = "pointillism"
    )
    markers <- as(import(file), "DataFrame")
    gene2symbol <- Gene2Symbol(seurat)
    keep <- markers[["geneID"]] %in% gene2symbol[["geneID"]]
    expect_true(any(keep))
    markers <- markers[keep, , drop = FALSE]
    x <- CellTypeMarkers(
        object = markers,
        gene2symbol = gene2symbol
    )
    expect_is(x, "CellTypeMarkers")
})
