context("counts")

test_that("Raw counts", {
    seurat <- as(sce, "Seurat")
    expect_identical(
        object = counts(seurat),
        expected = counts(sce)
    )

    cds <- as(sce, "cell_data_set")
    expect_identical(
        object = counts(cds),
        expected = counts(sce)
    )
})
