context("Coercion Methods")

test_that("Coerce SingleCellExperiment to seurat", {
    x <- as(sce_small, "seurat")
    expect_is(x, "seurat")
    # Check slotted count integrity
    counts <- counts(x)
    expect_is(counts, "dgCMatrix")
    expect_identical(
        dim(counts),
        dim(sce_small)
    )
})
