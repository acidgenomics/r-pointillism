context("Seurat")



test_that("SingleCellExperiment", {
    object <- Seurat::pbmc_small

    # assay
    expect_s4_class(assay(object), "dgCMatrix")

    # assayNames
    expect_identical(
        assayNames(object),
        c("counts", "logcounts")
    )

    # assays
    expect_s4_class(assays(object), "SimpleList")

    # colData
    expect_s4_class(colData(object), "DataFrame")

    # colData assignment
    x <- object
    colData(x)[["testthat"]] <- factor("XXX")
    expect_identical(
        levels(colData(x)[["testthat"]]),
        "XXX"
    )

    # colnames
    expect_is(colnames(object), "character")

    # counts
    expect_identical(counts(object), assay(object))

    # gene2symbol
    expect_warning(gene2symbol(object))
    expect_null(suppressWarnings(gene2symbol(object)))

    # interestingGroups
    expect_null(interestingGroups(object))

    # metadata
    expect_is(metadata(object), "list")

    # metadata assignment
    x <- object
    metadata(x)[["testthat"]] <- "XXX"
    expect_identical(
        metadata(x),
        list(testthat = "XXX")
    )

    # reducedDims
    x <- reducedDims(object)
    expect_s4_class(x, "SimpleList")
    expect_identical(names(x), c("PCA", "TSNE"))

    # rowData
    expect_s4_class(rowData(object), "DataFrame")

    # rownames
    expect_is(rownames(object), "character")

    # rowRanges
    expect_s4_class(rowRanges(object), "CompressedGRangesList")

    # sampleData
    expect_null(sampleData(object))

    # sampleNames
    expect_identical(sampleNames(object), character())
})



test_that("gene2symbol : seurat", {
    x <- gene2symbol(seurat_small)
    expect_is(x, "data.frame")
})

test_that("interestingGroups : seurat", {
    expect_null(interestingGroups(seurat_small))
})

test_that("interestingGroups<- : seurat", {
    # FIXME This is failing because of metadata assignment
    x <- seurat_small
    interestingGroups(x) <- "sampleName"
    expect_identical(
        interestingGroups(x),
        "sampleName"
    )
    x <- Seurat::pbmc_small
    expect_error(interestingGroups(x) <- "sampleName")
})



# metrics ======================================================================
test_that("metrics : seurat", {
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
            "orig.ident",
            "res.0.4",
            "res.0.8",
            "res.1.2",
            "sampleID",
            "sampleName"
        )
    )
})
