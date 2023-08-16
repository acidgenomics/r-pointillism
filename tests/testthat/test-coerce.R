test_that("Seurat to SingleCellExperiment", {
    object <- objs[["Seurat"]]
    x <- as(object, "SingleCellExperiment")
    expect_s4_class(x, "SingleCellExperiment")
    expect_identical(
        object = assayNames(x),
        expected = c("counts", "logcounts")
    )
    expect_identical(
        object = reducedDimNames(x),
        expected = c("PCA", "TSNE", "UMAP")
    )
    expect_identical(
        object = colnames(reducedDim(x, type = "UMAP")),
        expected = c("UMAP_1", "UMAP_2")
    )
    expect_identical(
        object = colnames(colData(x)),
        expected = c(
            "orig.ident",
            "nCount_RNA",
            "nFeature_RNA",
            "RNA_snn_res.0.8",
            "letter.idents",
            "groups",
            "RNA_snn_res.1",
            "ident"
        )
    )
    expect_identical(
        object = colnames(mcols(rowRanges(x))),
        expected = c(
            "broadClass",
            "entrezId",
            "geneBiotype",
            "geneId",
            "geneName",
            "seqCoordSystem"
        )
    )
    expect_named(
        object = metadata(x),
        expected = c(
            "scaleData",
            "variableFeatures"
        )
    )
})

test_that("SingleCellExperiment to Seurat", {
    object <- objs[["SingleCellExperiment"]]
    x <- as(object, "Seurat")
    expect_s4_class(x, "Seurat")
    counts <- counts(x)
    expect_s4_class(counts, "dgCMatrix")
    expect_identical(dim(counts), dim(object))
})

test_that("SCE-Seurat interconversion with subsetting", {
    a <- objs[["SingleCellExperiment"]]
    colDataNames <- colnames(colData(a))
    ## Coerce to Seurat.
    b <- as(a, "Seurat")
    expect_s4_class(b, "Seurat")
    ## FIXME This check is failing with Seurat 5.
    expect_identical(
        object = colnames(colData(b)),
        expected = colDataNames
    )
    ## Coerce back to SingleCellExperiment.
    c <- as(b, "SingleCellExperiment")
    expect_s4_class(c, "SingleCellExperiment")
    expect_identical(
        object = colnames(colData(c)),
        expected = colDataNames
    )
    ## Subset to contain n-1 genes, n-1 cells.
    nr <- nrow(c) - 1L
    nc <- ncol(c) - 1L
    d <- c[seq_len(nr), seq_len(nc)]
    expect_s4_class(d, "SingleCellExperiment")
    expect_identical(dim(d), c(nr, nc))
    ## Coerce back to Seurat.
    e <- as(d, "Seurat")
    expect_s4_class(e, "Seurat")
    expect_identical(dim(counts(e)), c(nr, nc))
    expect_identical(
        object = colnames(colData(e)),
        expected = colDataNames
    )
})
