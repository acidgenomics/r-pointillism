## NOTE Need to cover "seurat-clustering"
## NOTE Need to cover "seurat-markers"

render <- rmarkdown::render
seuratFile <- system.file(
    "data",
    "Seurat.rda",
    package = "AcidTest",
    mustWork = TRUE
)
tmpdir <- AcidBase::tempdir2()

test_that("Seurat per cluster analysis", {
    skeleton <- system.file(
        "rmarkdown",
        "templates",
        "seurat-per-cluster-analysis",
        "skeleton",
        "skeleton.Rmd",
        package = .pkgName,
        mustWork = TRUE
    )
    input <- tempfile(
        pattern = "render",
        tmpdir = tmpdir,
        fileext = ".Rmd"
    )
    file.copy(from = skeleton, to = input, overwrite = TRUE)
    out <- render(
        input = input,
        output_format = "md_document",
        clean = TRUE,
        params = list("seurat_file" = seuratFile),
        quiet = TRUE
    )
    expect_true(file.exists(out))
})

AcidBase::unlink2(tmpdir)
