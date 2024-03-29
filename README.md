# pointillism

![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)

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
        "https://r.acidgenomics.com",
        BiocManager::repositories()
    ),
    dependencies = TRUE
)
```

## Supported data classes

pointillism currently supports these S4 single-cell container classes:

- [SingleCellExperiment][]
- [Seurat][]
- [monocle3][] `cell_data_set`

## Markers

Cell-cycle and cell-type markers are stored internally inside the package.
Refer to `inst/extdata/` for the source CSV files.

## Troubleshooting

### Maximal number of DLLs reached

```txt
Error: package or namespace load failed for 'pointillism'
in dyn.load(file, DLLpath = DLLpath, ...):
  maximal number of DLLs reached...
```

Depending on your operating system, you may encounter this error about hitting
the DLL limit in [R][]. This issue is becoming more common as RNA-seq analysis
packages grow increasingly complex. Luckily, we can configure [R][] to increase
the DLL limit. Append this line to your `~/.Renviron` file:

```sh
R_MAX_NUM_DLLS=150
```

For more information on this issue, consult `help("dyn.load")` in the [R][]
documentation. The number of loaded DLLs in an [R][] session can be obtained
with `getLoadedDLLs()`.

## References

The papers and software cited in our workflows are available as a
[shared library](https://paperpile.com/shared/5PLRi1) on [Paperpile][].

[monocle3]: https://cole-trapnell-lab.github.io/monocle3/
[paperpile]: https://paperpile.com/
[r]: https://www.r-project.org/
[seurat]: https://satijalab.org/seurat/
[singlecellexperiment]: https://bioconductor.org/packages/SingleCellExperiment/
