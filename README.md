# pointillism

[![Travis CI](https://travis-ci.org/steinbaugh/pointillism.svg?branch=master)](https://travis-ci.org/steinbaugh/pointillism)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

[R][] package for for single-cell RNA-seq clustering analysis.

## Installation

This is an [R][] package.

### [Bioconductor][]

We recommend installing the package with [BiocManager][].

```r
if (!require("BiocManager")) {
    install.packages("BiocManager")
}
BiocManager::install("remotes")
BiocManager::install("acidgenomics/pointillism")
```

## Supported data classes

[pointillism][] currently supports these S4 single-cell container classes:

- [SingleCellExperiment][]
- [Seurat][] (v3)

**monocle v3**: Support for the [monocle][] `CellDataSet` class will be added when [monocle][] v3 becomes available on [Bioconductor][].

## Markers

Shared [cell-cycle markers][] and [cell-type markers][] are available on [Google Sheets][]. Contact [Michael Steinbaugh][] if you'd like to contribute to this list.

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

[Bioconductor]: https://bioconductor.org/
[BiocManager]: https://cran.r-project.org/package=BiocManager
[Cell-cycle markers]: https://docs.google.com/spreadsheets/d/1qA5ktYeimNGpZF1UPSQZATbpzEqgyxN6daoMOjv6YYw
[Cell-type markers]: https://docs.google.com/spreadsheets/d/1vGNU2CCxpaoTCLvzOxK1hf5gjULrf2-CpgCp9bOfGJ0
[conda]: https://conda.io/
[CRAN]: https://cran.r-project.org/  "The Comprehensive R Archive Network"
[Google Sheets]: https://www.google.com/sheets
[Michael Steinbaugh]: https://mike.steinbaugh.com/
[monocle]: http://cole-trapnell-lab.github.io/monocle-release/
[Paperpile]: https://paperpile.com/
[pointillism]: https://pointillism.acidgenomics.com/
[R]: https://www.r-project.org/
[Seurat]: https://satijalab.org/seurat/
[SingleCellExperiment]: https://bioconductor.org/packages/SingleCellExperiment/
