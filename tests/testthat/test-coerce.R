context("coerce")

test_that("Seurat to SingleCellExperiment", {
    x <- as(seurat, "SingleCellExperiment")
    expect_s4_class(x, "SingleCellExperiment")
    expect_identical(
        object = assayNames(x),
        expected = c("counts", "logcounts")
    )
    expect_identical(
        object = reducedDimNames(x),
        expected = c("pca", "tsne", "umap")
    )
    expect_identical(
        object = names(metadata(x)),
        expected = c(
            "scaleData",
            "variableFeatures"
        )
    )
})

test_that("SingleCellExperiment to Seurat", {
    x <- as(sce, "Seurat")
    expect_is(x, "Seurat")
    ## Check slotted count integrity.
    counts <- counts(x)
    expect_is(counts, "dgCMatrix")
    expect_identical(dim(counts), dim(sce))
})

## FIXME NEED TO TEST FOR MESSED UP COLUMN NAMES HERE...
## CANT CAMELCASE THE SCE, OTHERWISE WE GET DUPES IN SEURAT ARGH...
test_that("SCE-seurat interconversion with subsetting", {
    a <- sce
    ## Coerce to seurat.
    b <- as(a, "Seurat")
    expect_s4_class(b, "Seurat")
    ## Coerce back to SCE.
    c <- as(b, "SingleCellExperiment")
    expect_s4_class(c, "SingleCellExperiment")
    ## Subset to contain n-1 genes, n-1 cells.
    nr <- nrow(c) - 1L
    nc <- ncol(c) - 1L
    d <- c[seq_len(nr), seq_len(nc)]
    expect_s4_class(d, "SingleCellExperiment")
    expect_identical(dim(d), c(nr, nc))
    ## Coerce back to seurat.
    e <- as(d, "Seurat")
    expect_s4_class(e, "Seurat")
    expect_identical(dim(counts(e)), c(nr, nc))
})
