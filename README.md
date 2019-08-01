# pointillism

[![Repo status: active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Travis CI build status](https://travis-ci.com/acidgenomics/pointillism.svg?branch=master)](https://travis-ci.com/acidgenomics/pointillism)

[R][] package for for single-cell RNA-seq clustering analysis.

## Installation

### [R][] method

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
Sys.setenv(R_REMOTES_UPGRADE = "always")
# Set `GITHUB_PAT` in `~/.Renviron` if you get a rate limit error.
remotes::install_github("acidgenomics/pointillism")
```

Here's how to update to the latest version on GitHub:

```r
Sys.setenv(R_REMOTES_UPGRADE = "always")
remotes::update_packages()
```

Always check that your Bioconductor installation is valid before proceeding.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::valid()
```

## Supported data classes

[pointillism][] currently supports these S4 single-cell container classes:

- [SingleCellExperiment][]
- [Seurat][] (v3)

We're working on adding support for monocle3 `cell_data_set` class for all plotting functions in the current v0.4 release series.

## Markers

Cell-cycle and cell-type markers are stored internally inside the package. Refer to `inst/extdata/` for the source CSV files.

## Troubleshooting

### Maximal number of DLLs reached

```
Error: package or namespace load failed for 'pointillism' in dyn.load(file, DLLpath = DLLpath, ...):
  maximal number of DLLs reached...
```

Depending on your operating system, you may encounter this error about hitting the DLL limit in [R][]. This issue is becoming more common as RNA-seq analysis packages grow increasingly complex. Luckily, we can configure [R][] to increase the DLL limit. Append this line to your `~/.Renviron` file:

```
R_MAX_NUM_DLLS=150
```

For more information on this issue, consult `help("dyn.load")` in the [R][] documentation. The number of loaded DLLs in an [R][] session can be obtained with `getLoadedDLLs()`.

## References

The papers and software cited in our workflows are available as a [shared library](https://paperpile.com/shared/5PLRi1) on [Paperpile][].

[BiocManager]: https://cran.r-project.org/package=BiocManager
[Bioconductor]: https://bioconductor.org/
[CRAN]: https://cran.r-project.org/  "The Comprehensive R Archive Network"
[Michael Steinbaugh]: https://mike.steinbaugh.com/
[Paperpile]: https://paperpile.com/
[R]: https://www.r-project.org/
[Seurat]: https://satijalab.org/seurat/
[SingleCellExperiment]: https://bioconductor.org/packages/SingleCellExperiment/
[conda]: https://conda.io/
[monocle]: http://cole-trapnell-lab.github.io/monocle-release/
[pointillism]: https://pointillism.acidgenomics.com/
