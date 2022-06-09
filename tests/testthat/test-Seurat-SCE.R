object <- objs[["Seurat"]]

test_that("Gene2Symbol", {
    x <- Gene2Symbol(object)
    expect_is(x, "Gene2Symbol")
})

test_that("assay", {
    expect_s4_class(assay(object), "sparseMatrix")
})

test_that("assayNames", {
    expect_identical(
        object = assayNames(object),
        expected = c("counts", "logcounts")
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

test_that("interestingGroups", {
    expect_null(interestingGroups(object))
})

test_that("interestingGroups<-", {
    expect_silent(
        interestingGroups(object) <- "sampleName"
    )
    expect_error(
        interestingGroups(object) <- "orig.ident"
    )
    expect_error(
        interestingGroups(object) <- "XXX"
    )
    interestingGroups(object) <- "sampleName"
    expect_identical(
        object = interestingGroups(object),
        expected = "sampleName"
    )
})

test_that("metadata", {
    expect_is(metadata(object), "list")
    ## Assignment method.
    metadata(object)[["testthat"]] <- "XXX"
    expect_identical(
        object = metadata(object)[["testthat"]],
        expected = "XXX"
    )
})

test_that("metrics", {
    x <- metrics(object)
    expect_is(x, "DataFrame")
})

test_that("reducedDims", {
    x <- reducedDims(object)
    expect_s4_class(x, "SimpleList")
    expect_named(x, c("PCA", "TSNE", "UMAP"))
})

test_that("rowData", {
    expect_s4_class(rowData(object), "DataFrame")
})

test_that("rownames", {
    expect_type(rownames(object), "character")
})

test_that("rowRanges", {
    expect_s4_class(rowRanges(object), "GenomicRanges")
})

test_that("sampleData", {
    expect_identical(
        object = sampleData(object),
        expected = DataFrame(
            "sampleName" = factor("unknown"),
            "interestingGroups" = factor("unknown"),
            row.names = "unknown"
        )
    )
})

test_that("sampleNames", {
    expect_identical(
        object = sampleNames(object),
        expected = c(unknown = "unknown")
    )
})
