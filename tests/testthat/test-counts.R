test_that("Seurat", {
    expect_identical(
        object = counts(objs[["Seurat"]]),
        expected = counts(objs[["SingleCellExperiment"]])
    )
})
