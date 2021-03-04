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
            "seqCoordSystem",
            "vstMean",
            "vstVariance",
            "vstVarianceExpected",
            "vstVarianceStandardized",
            "vstVariable"
        )
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
    colDataNames <- colnames(colData(a))
    ## Coerce to seurat.
    b <- as(a, "Seurat")
    expect_s4_class(b, "Seurat")
    expect_identical(
        object = colnames(colData(b)),
        expected = colDataNames
    )
    ## Coerce back to SCE.
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
    ## Coerce back to seurat.
    e <- as(d, "Seurat")
    expect_s4_class(e, "Seurat")
    expect_identical(dim(counts(e)), c(nr, nc))
    expect_identical(
        object = colnames(colData(e)),
        expected = colDataNames
    )
})
