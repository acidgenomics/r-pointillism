context("Seurat as SingleCellExperiment")

pbmc_small <- Seurat::pbmc_small



test_that("assay", {
    expect_s4_class(assay(pbmc_small), "dgCMatrix")
})



test_that("assayNames", {
    expect_identical(
        assayNames(pbmc_small),
        c("counts", "logcounts")
    )
})



test_that("assays", {
    expect_s4_class(assays(pbmc_small), "SimpleList")
})



test_that("colData", {
    expect_s4_class(colData(pbmc_small), "DataFrame")
})



test_that("colData<-", {
    x <- pbmc_small
    colData(x)[["testthat"]] <- factor("XXX")
    expect_identical(
        levels(colData(x)[["testthat"]]),
        "XXX"
    )
})



test_that("colnames", {
    expect_is(colnames(pbmc_small), "character")
})



test_that("counts", {
    expect_identical(counts(pbmc_small), assay(pbmc_small))
})



test_that("gene2symbol", {
    expect_warning(Gene2Symbol(pbmc_small))
    expect_null(suppressWarnings(Gene2Symbol(pbmc_small)))

    x <- Gene2Symbol(seurat_small)
    expect_is(x, "data.frame")
})



test_that("interestingGroups", {
    expect_null(interestingGroups(pbmc_small))
})



test_that("interestingGroups<-", {
    x <- pbmc_small
    # We're requiring `sampleData()` return here, which requires `sampleID`
    # and `sampleName` columns in `colData()`.
    expect_error(
        interestingGroups(x) <- "sampleName",
        "colData"
    )
    interestingGroups(x) <- "orig.ident"
    expect_identical(
        interestingGroups(x),
        "orig.ident"
    )

    x <- seurat_small
    interestingGroups(x) <- "sampleName"
    expect_identical(
        interestingGroups(x),
        "sampleName"
    )
    expect_error(
        interestingGroups(x) <- "XXX",
        "colData"
    )
})



test_that("metadata", {
    expect_is(metadata(pbmc_small), "list")

    # metadata assignment
    x <- pbmc_small
    metadata(x)[["testthat"]] <- "XXX"
    expect_identical(
        metadata(x),
        list(testthat = "XXX")
    )
})



test_that("metrics", {
    # Check for camel case sanitization in `as()` coercion method.
    expect_identical(
        sort(colnames(metrics(seurat_small))),
        c(
            "batch",
            "expLibSize",
            "group",
            "ident",
            "interestingGroups",
            "log10GenesPerUMI",
            "mitoRatio",
            "nCoding",
            "nGene",
            "nMito",
            "nUMI",
            "origIdent",
            "res0.4",
            "res0.8",
            "res1.2",
            "sampleID",
            "sampleName"
        )
    )
})



test_that("reducedDims", {
    x <- reducedDims(pbmc_small)
    expect_s4_class(x, "SimpleList")
    expect_identical(names(x), c("PCA", "TSNE"))
})



test_that("rowData", {
    expect_s4_class(rowData(pbmc_small), "DataFrame")
})



test_that("rownames", {
    expect_is(rownames(pbmc_small), "character")
})



test_that("rowRanges", {
    expect_s4_class(rowRanges(pbmc_small), "CompressedGRangesList")
})



test_that("sampleData", {
    expect_error(sampleData(pbmc_small), "sampleData")

    x <- sampleData(sce_small)
    expect_identical(
        rownames(x),
        c("group1", "group2")
    )
    expect_identical(
        colnames(x),
        c(
            "batch",
            "group",
            "sampleID",
            "sampleName",
            "interestingGroups"
        )
    )
})



test_that("sampleNames", {
    expect_error(
        sampleNames(pbmc_small),
        "sampleData"
    )
    expect_identical(
        sampleNames(seurat_small),
        c(
            group1 = "group1",
            group2 = "group2"
        )
    )
})
