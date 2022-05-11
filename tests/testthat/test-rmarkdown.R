## NOTE Need to cover "seurat-clustering"
## NOTE Need to cover "seurat-markers"



context("R Markdown templates")

render <- rmarkdown::render

seuratFile <- system.file(
    "data",
    "seurat.rda",
    package = "pointillism",
    mustWork = TRUE
)

tmpdir <- file.path(tempdir(), "rmarkdown-render")
unlink(tmpdir, recursive = TRUE)
dir.create(tmpdir, recursive = TRUE)

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

unlink(tmpdir, recursive = TRUE)
