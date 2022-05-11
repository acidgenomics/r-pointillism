context("R Markdown")

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
    ## FIXME This is hanging at the end, need to figure out why.
    ## FIXME Consider changing output format back to HTML.
    ## FIXME pandoc is going crazy on this example, may be too memory intensive...
    out <- rmarkdown::render(
        input = input,
        output_format = "md_document",
        clean = TRUE,
        params = list(
            "seurat_file" = system.file(
                "data",
                "Seurat.rda",
                package = "pointillism",
                mustWork = TRUE
            )
        ),
        quiet = TRUE
    )
    expect_true(file.exists(out))
})

unlink(tmpdir, recursive = TRUE)
