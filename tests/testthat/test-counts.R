context("counts")

test_that("Seurat", {
    seurat <- as(sce, "Seurat")
    expect_identical(
        object = counts(seurat),
        expected = counts(sce)
    )
})

## test_that("cell_data_set", {
##     cds <- as(sce, "cell_data_set")
##     expect_identical(
##         object = counts(cds),
##         expected = counts(sce)
##     )
## })
