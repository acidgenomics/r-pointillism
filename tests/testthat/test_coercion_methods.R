context("Coercion Methods")

test_that("SingleCellExperiment to seurat", {
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

test_that("SCE/seurat interconversion with subsetting", {
    a <- sce_small
    # Coerce to seurat
    b <- as(a, "seurat")
    expect_s4_class(b, "seurat")
    # Coerce back to SCE
    c <- as(b, "SingleCellExperiment")
    expect_s4_class(c, "SingleCellExperiment")
    # Subset to contain 100 genes, 10 cells
    d <- c[seq_len(100L), seq_len(10L)]
    expect_s4_class(d, "SingleCellExperiment")
    expect_identical(dim(d), c(100L, 10L))
    # Coerce back to seurat
    e <- as(d, "seurat")
    expect_s4_class(e, "seurat")
    expect_identical(dim(counts(e)), c(100L, 10L))
})
