context("counts")

test_that("Seurat", {
    expect_identical(
        object = counts(Seurat),
        expected = counts(SingleCellExperiment)
    )
})

## > test_that("cell_data_set", {
## >     expect_identical(
## >         object = counts(cell_data_set),
## >         expected = counts(SingleCellExperiment)
## >     )
## > })
