test_that("Seurat", {
    expect_identical(
        object = counts(objs[["Seurat"]]),
        expected = counts(objs[["SingleCellExperiment"]])
    )
})

## > test_that("cell_data_set", {
## >     expect_identical(
## >         object = counts(objs[["cell_data_set"]]),
## >         expected = counts(objs[["SingleCellExperiment"]])
## >     )
## > })
