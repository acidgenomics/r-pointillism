# pointillism

[![Repo status: active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Travis CI build status](https://travis-ci.com/acidgenomics/pointillism.svg?branch=master)](https://travis-ci.com/acidgenomics/pointillism)

[R][] package for for single-cell RNA-seq clustering analysis.

## Installation

### [R][] method

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
install.packages(
    pkgs = "pointillism",
    repos = c(
        "r.acidgenomics.com",
        BiocManager::repositories()
    )
)
```

## Supported data classes

[pointillism][] currently supports these S4 single-cell container classes:

- [SingleCellExperiment][]
- [Seurat][]
- [monocle3][] `cell_data_set`

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

[biocmanager]: https://cran.r-project.org/package=BiocManager
[bioconductor]: https://bioconductor.org/
[conda]: https://conda.io/
[cran]: https://cran.r-project.org/  "The Comprehensive R Archive Network"
[michael steinbaugh]: https://mike.steinbaugh.com/
[monocle3]: https://cole-trapnell-lab.github.io/monocle3/
[paperpile]: https://paperpile.com/
[pointillism]: https://pointillism.acidgenomics.com/
[r]: https://www.r-project.org/
[seurat]: https://satijalab.org/seurat/
[singlecellexperiment]: https://bioconductor.org/packages/SingleCellExperiment/
