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

    # sampleNames
    expect_identical(sampleNames(object), character())
})



test_that("gene2symbol : seurat", {
    x <- gene2symbol(seurat_small)
    expect_is(x, "data.frame")
    expect_identical(colnames(x), colnames)
})

test_that("interestingGroups : seurat", {
    expect_identical(
        interestingGroups(seurat_small),
        "sampleName"
    )
    expect_identical(
        interestingGroups(seurat_small),
        "sampleName"
    )
})

test_that("interestingGroups<- : seurat", {
    interestingGroups(indrops_small) <- "sampleName"
    expect_identical(
        interestingGroups(indrops_small),
        "sampleName"
    )
    interestingGroups(seurat_small) <- "sampleName"
    expect_identical(
        interestingGroups(seurat_small),
        "sampleName"
    )
    x <- Seurat::pbmc_small
    expect_error(interestingGroups(Seurat::pbmc_small) <- "sampleName")
})



# metrics ======================================================================
test_that("metrics : seurat", {
    # Check that metrics accessor data matches meta.data slot
    x <- metrics(seurat_small)
    y <- seurat_small@meta.data
    x <- x[, colnames(y)]
    expect_identical(x, y)
})



test_that("sampleData : seurat", {
    # Return all columns
    x <- sampleData(seurat_small)
    expect_identical(
        lapply(x, class),
        list(
            sampleName = "factor",
            description = "factor",
            index = "factor",
            interestingGroups = "factor"
        )
    )

    # Return NULL for other seurat objects
    expect_identical(sampleData(Seurat::pbmc_small), NULL)
})
