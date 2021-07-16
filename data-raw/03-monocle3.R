## monocle3 cell_data_set example.
## Updated 2021-07-16.

## This is currently failing to save with Bioconductor 3.10, due to changes in
## SingleCellExperiment that now cause validity checks to fail.
## See related issue on GitHub:
## https://github.com/cole-trapnell-lab/monocle3/issues/246

## See also:
## - https://cole-trapnell-lab.github.io/monocle3/monocle3_docs/
## - https://github.com/cole-trapnell-lab/monocle3/blob/master/examples/
##       c_elegans_L2.R
## - https://github.com/cole-trapnell-lab/monocle3/blob/master/examples/
##       c_elegans_embryo.R

suppressPackageStartupMessages({
    library(usethis)     # 2.0.1
    library(magrittr)    # 1.5
    library(lobstr)      # 1.1.1
    library(reticulate)  # 1.20
    library(goalie)      # 0.5.1
    library(basejump)    # 0.14.19
    library(monocle3)    # 1.0.0
})

cores <- getOption("mc.cores")

virtualenv_list()
## [1] "r-reticulate"

use_virtualenv(virtualenv = "r-reticulate", required = TRUE)

py_config()
## python:         /opt/koopa/opt/virtualenvs/r-reticulate/bin/python
## libpython:      /Library/Frameworks/Python.framework/Versions/3.9/lib/python3.9/config-3.9-darwin/libpython3.9.dylib
## pythonhome:     /opt/koopa/opt/virtualenvs/r-reticulate:/opt/koopa/opt/virtualenvs/r-reticulate
## version:        3.9.6 (v3.9.6:db3ff76da1, Jun 28 2021, 11:14:58)  [Clang 12.0.5 (clang-1205.0.22.9)]
## numpy:          /opt/koopa/opt/virtualenvs/r-reticulate/lib/python3.9/site-packages/numpy
## numpy_version:  1.21.0

## Restrict object size to 2 MB.
## Use `pryr::object_size()` instead of `utils::object.size()`.
limit <- structure(2e6, class = "object_size")

## Load data ===================================================================
prefix <- "http://staff.washington.edu/hpliner/data"

## Don't time out after 60 seconds here.
## The University of Washington staff website is currently pretty slow.
options("timeout" = 600L)
expression_data <- import(file.path(prefix, "cao_l2_expression.rds"))
## > cell_metadata <- import(file.path(prefix, "cao_l2_colData.rds"))
gene_metadata <- import(file.path(prefix, "cao_l2_rowData.rds"))

head(colnames(expression_data))
## Don't allow "-" in the names, not syntactically valid.
## [1] "cele-001-001.CATGACTCAA" "cele-001-001.AAGACGGCCA"
## [3] "cele-001-001.GCCAACGCCA" "cele-001-001.ATAGGAGTAC"
## [5] "cele-001-001.CTCGTCTAGG" "cele-001-001.AAGTTGCCAT"

## Make the dimnames valid before proceeding.
expression_data <- makeDimnames(expression_data)

assert(
    hasValidDimnames(expression_data),
    hasValidDimnames(gene_metadata)
)

## If you hit this error below, reinstall monocle3.
## > remotes::install_github("cole-trapnell-lab/monocle3", force = TRUE)
##
## Error in checkSlotAssignment(object, name, value) :
## 'reducedDims' is not a slot in class "SingleCellExperiment"
## Calls: new_cell_data_set ... .valid.Vector.length -> as -> asMethod ->
##     slot<- -> checkSlotAssignment
##
## https://github.com/cole-trapnell-lab/monocle3/issues/246

## Example cell metadata contains cluster mappings, so skip loading that.
cds <- new_cell_data_set(
    expression_data = expression_data,
    ## cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
)

## Subset to include only the top cells and genes by number of reads.
counts <- counts(cds)

topGenes <- counts %>%
    rowSums() %>%
    sort(decreasing = TRUE) %>%
    head(n = 200L) %>%
    names()
topCells <- counts %>%
    colSums() %>%
    sort(decreasing = TRUE) %>%
    head(n = 200L) %>%
    names()

cds <- cds[topGenes, topCells]

## Relevel the factors to decrease object size.
rowData(cds) <- droplevels(rowData(cds))
colData(cds) <- droplevels(colData(cds))

slotNames(cds)
##  [1] "preprocess_aux"      "reduce_dim_aux"      "principal_graph_aux"
##  [4] "principal_graph"     "clusters"            "int_elementMetadata"
##  [7] "int_colData"         "int_metadata"        "rowRanges"
## [10] "colData"             "assays"              "NAMES"
## [13] "elementMetadata"     "metadata"

colnames(colData(cds))
## [1] "cell"        "Size_Factor"

## Pre-process the data.
## Note that log normalization of data has come under question as appropriate.
cds <- preprocess_cds(cds)

## Reduce dimensionality.
## Note: reduce_dimension will produce slightly different output each time you
## run it unless you set 'umap.fast_sgd = FALSE' and 'cores = 1'
cds <- reduce_dimension(
    cds = cds,
    reduction_method = "UMAP",
    preprocess_method = "PCA",
    umap.fast_sgd = TRUE,
    cores = cores
)

## Group cells into clusters.
cds <- cluster_cells(
    cds = cds,
    reduction_method = "UMAP",
    k = 20L,
    louvain_iter = 1L,
    partition_qval = 0.05,
    weight = FALSE,
    resolution = c(10L ^ seq(from = -6L, to = -1L)),
    verbose = TRUE
)

## Clusters are currently located in cds@clusters
## cds@clusters[["UMAP"]]

## Order cells in pseudotime along a trajectory.
## This step takes a while for large datasets.
cds <- learn_graph(
    cds = cds,
    use_partition = TRUE,
    close_loop = TRUE,
    verbose = TRUE
)

stopifnot(obj_size(cds) < limit)

cell_data_set <- cds
use_data(cell_data_set, compress = "xz", overwrite = TRUE)
