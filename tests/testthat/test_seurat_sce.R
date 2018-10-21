context("Seurat as SingleCellExperiment")

data(seurat_small)
object <- seurat_small



test_that("assay", {
    expect_s4_class(assay(object), "sparseMatrix")
})



test_that("assayNames", {
    expect_identical(
        assayNames(object),
        c("counts", "logcounts")
    )
})



test_that("assays", {
    expect_s4_class(assays(object), "SimpleList")
})



test_that("colData", {
    expect_s4_class(colData(object), "DataFrame")
})



test_that("colData<-", {
    x <- object
    colData(x)[["testthat"]] <- factor("XXX")
    expect_identical(
        levels(colData(x)[["testthat"]]),
        "XXX"
    )
})



test_that("colnames", {
    expect_is(colnames(object), "character")
})



test_that("counts", {
    expect_identical(counts(object), assay(object))
})



test_that("Gene2Symbol", {
    x <- Gene2Symbol(seurat_small)
    expect_is(x, "Gene2Symbol")
})



test_that("interestingGroups", {
    expect_null(interestingGroups(object))
})



test_that("interestingGroups<-", {
    expect_error(
        object = interestingGroups(object) <- "orig.ident",
        regexp = "sampleData"
    )
    expect_error(
        object = interestingGroups(object) <- "XXX",
        regexp = "sampleData"
    )

    expect_silent(
        object = interestingGroups(object) <- "sampleName"
    )
    interestingGroups(object) <- "sampleName"
    expect_identical(
        object = interestingGroups(object),
        expected = "sampleName"
    )
})



test_that("metadata", {
    expect_is(metadata(object), "list")
    # Assignment method.
    metadata(object)[["testthat"]] <- "XXX"
    expect_identical(
        object = metadata(object)[["testthat"]],
        expected = "XXX"
    )
})



test_that("metrics", {
    expect_identical(
        sort(colnames(metrics(object))),
        c(
            "cellID",
            "ident",
            "interestingGroups",
            "nGene",
            "nUMI",
            "orig.ident",
            "res.0.8",
            "res.1",
            "sampleID",
            "sampleName"
        )
    )
})



test_that("reducedDims", {
    x <- reducedDims(object)
    expect_s4_class(x, "SimpleList")
    expect_identical(names(x), c("PCA", "TSNE", "UMAP"))
})



test_that("rowData", {
    expect_s4_class(rowData(object), "DataFrame")
})



test_that("rownames", {
    expect_type(rownames(object), "character")
})



test_that("rowRanges", {
    expect_s4_class(rowRanges(object), "GRanges")
})



# FIXME Interesting groups column is named in `sampleData` return.
test_that("sampleData", {
    expect_identical(
        object = sampleData(object),
        expected = DataFrame(
            sampleName = factor("unknown"),
            interestingGroups = factor("unknown"),
            row.names = "unknown"
        )
    )
})



test_that("sampleNames", {
    expect_error(
        sampleNames(object),
        "sampleData"
    )
    expect_identical(
        sampleNames(object),
        c(
            group1 = "group1",
            group2 = "group2"
        )
    )
})
