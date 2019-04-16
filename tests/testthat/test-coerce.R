# FIXME Need to test seurat to SCE.
# FIXME Need to test Markers to tbl_df.

context("coerce")

test_that("Seurat to SingleCellExperiment", {
    x <- as(seurat
})

test_that("SingleCellExperiment to Seurat", {
    x <- as(sce, "seurat")
    expect_is(x, "seurat")
    # Check slotted count integrity.
    counts <- counts(x)
    expect_is(counts, "dgCMatrix")
    expect_identical(dim(counts), dim(sce))
})

test_that("SCE-seurat interconversion with subsetting", {
    a <- sce

    # Coerce to seurat.
    b <- as(a, "seurat")
    expect_s4_class(b, "seurat")

    # Coerce back to SCE.
    c <- as(b, "SingleCellExperiment")
    expect_s4_class(c, "SingleCellExperiment")

    # Subset to contain n-1 genes, n-1 cells.
    nr <- nrow(c) - 1L
    nc <- ncol(c) - 1L
    d <- c[seq_len(nr), seq_len(nc)]
    expect_s4_class(d, "SingleCellExperiment")
    expect_identical(dim(d), c(nr, nc))

    # Coerce back to seurat.
    e <- as(d, "seurat")
    expect_s4_class(e, "seurat")
    expect_identical(dim(counts(e)), c(nr, nc))
})
