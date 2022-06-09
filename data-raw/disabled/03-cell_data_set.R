## monocle3 cell_data_set example.
##
## @note Updated 2022-06-09.
##
## @seealso
## - https://cole-trapnell-lab.github.io/monocle3/monocle3_docs/
## - https://github.com/cole-trapnell-lab/monocle3/blob/master/examples/
## c_elegans_L2.R
## - https://github.com/cole-trapnell-lab/monocle3/blob/master/examples/
## c_elegans_embryo.R
## - https://github.com/cole-trapnell-lab/monocle3/issues/246
## nolint start
suppressPackageStartupMessages({
    library(usethis)
    library(reticulate)
    library(goalie)
    library(basejump)
    library(monocle3)
    library(pointillism)
})
## nolint end
reticulate::use_condaenv(
    condaenv = "umap-learn@0.5.2",
    required = TRUE
)
limit <- structure(2e6L, class = "object_size")
prefix <- "http://staff.washington.edu/hpliner/data"
## Don't time out after 60 seconds here.
## The University of Washington staff website is currently pretty slow.
# > options("timeout" = 600L)
expressionData <-
    file.path(prefix, "cao_l2_expression.rds") |>
    import() |>
    makeDimnames()
cellMetadata <-
    file.path(prefix, "cao_l2_colData.rds") |>
    import() |>
    makeDimnames()
geneMetadata <-
    file.path(prefix, "cao_l2_rowData.rds") |>
    import() |>
    makeDimnames()
cds <- new_cell_data_set(
    expression_data = expressionData,
    cell_metadata = cellMetadata,
    gene_metadata = geneMetadata
)
## Subset to include only the top cells and genes by number of reads.
topGenes <-
    counts(cds) |>
    rowSums() |>
    sort(decreasing = TRUE) |>
    head(n = 200L) |>
    names()
topCells <-
    counts(cds) |>
    colSums() |>
    sort(decreasing = TRUE) |>
    head(n = 200L) |>
    names()
cds <- cds[topGenes, topCells]
## Relevel the factors to decrease object size.
rowData(cds) <- droplevels(rowData(cds))
colData(cds) <- droplevels(colData(cds))
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
    cores = getOption(x = "mc.cores")
)
## Group cells into clusters.
cds <- cluster_cells(
    cds = cds,
    reduction_method = "UMAP",
    k = 20L,
    louvain_iter = 1L,
    partition_qval = 0.05,
    weight = FALSE,
    resolution = c(10L^seq(from = -6L, to = -1L)),
    verbose = TRUE
)
## Clusters are currently located in cds@clusters
## > cds@clusters[["UMAP"]]
## Order cells in pseudotime along a trajectory.
## This step takes a while for large datasets.
cds <- learn_graph(
    cds = cds,
    use_partition = TRUE,
    close_loop = TRUE,
    verbose = TRUE
)
stopifnot(object.size(cds) < limit)
cell_data_set <- cds # nolint
use_data(cell_data_set, compress = "xz", overwrite = TRUE)
