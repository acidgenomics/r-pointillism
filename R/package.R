#' pointillism
#'
#' R package for for single-cell RNA-seq clustering analysis.
#'
#' @keywords internal
#'
#' @importClassesFrom Seurat seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom assertive.base assert_are_identical assert_all_are_not_na
#' @importFrom assertive.numbers assert_all_are_in_left_open_range
#'   assert_all_are_in_open_range assert_all_are_non_negative
#'   assert_all_are_positive
#' @importFrom assertive.properties assert_has_no_duplicates assert_has_names
#'   assert_has_rows assert_is_non_empty assert_is_of_length assert_is_scalar
#' @importFrom assertive.sets assert_are_disjoint_sets
#'   assert_are_intersecting_sets assert_is_subset
#' @importFrom assertive.types assert_is_a_bool assert_is_a_number
#'   assert_is_a_string assert_is_all_of assert_is_an_integer assert_is_any_of
#'   assert_is_character assert_is_data.frame assert_is_factor
#'   assert_is_function assert_is_list assert_is_matrix assert_is_numeric
#'   assert_is_tbl_df is_a_string
#' @importFrom assertthat assert_that validate_that
#' @importFrom basejump camel coerce
#'   convertGenesToSymbols gene2symbol import mapGenesToRownames
#'   mapGenesToSymbols markdownHeader matchArgsToDoCall matchInterestingGroups
#'   printString snake theme_midnight upperCamel
#' @importFrom BiocGenerics cbind counts counts<- do.call lapply sapply
#' @importFrom BiocParallel bpparam bpprogressbar bpprogressbar<- MulticoreParam
#'   SerialParam
#' @importFrom cowplot plot_grid
#' @importFrom dplyr arrange desc everything filter group_by group_vars
#'   left_join mutate mutate_at mutate_if n pull rename select slice summarize
#'   ungroup vars
#' @importFrom DESeq2 DESeqDataSet DESeq results
#' @importFrom edgeR calcNormFactors DGEList estimateDisp glmFit
#' @importFrom ggplot2 aes element_rect expand_limits facet_grid facet_wrap
#'   geom_hline geom_jitter geom_line geom_point geom_text geom_violin
#'   geom_vline ggplot guide_colorbar guides labs scale_color_gradient
#'   scale_color_viridis_c scale_radius scale_x_continuous scale_y_continuous
#'   theme
#' @importFrom goalie assertHasRownames assertIsAlpha assertIsAStringOrNULL
#'   assertIsAnImplicitInteger assertIsColorScaleContinuousOrNULL
#'   assertIsColorScaleDiscreteOrNULL assertIsFillScaleDiscreteOrNULL
#'   assertIsHeaderLevel assertIsImplicitInteger
#' @importFrom googlesheets gs_key gs_read
#' @importFrom magrittr %>% set_names set_rownames
#' @importFrom Matrix rowMeans rowSums
#' @importFrom methods as getMethod is new setAs show slot slot<- validObject
#' @importFrom pbapply pblapply
#' @importFrom rlang !! !!! := has_length sym syms
#' @importFrom S4Vectors as.data.frame as.matrix complete.cases head mcols
#'   mcols<- metadata metadata<- na.omit split tail
#' @importFrom scales percent pretty_breaks
#' @importFrom sessioninfo session_info
#' @importFrom Seurat Convert CreateSeuratObject
#' @importFrom SingleCellExperiment reducedDimNames reducedDims weights
#'   weights<-
#' @importFrom stats median model.matrix relevel
#' @importFrom SummarizedExperiment assay assayNames assays assays<- colData
#'   rowData rowRanges rowRanges<-
#' @importFrom tibble as_tibble column_to_rownames remove_rownames tibble
#' @importFrom tidyr gather
#' @importFrom utils capture.output globalVariables packageVersion
#' @importFrom zinbwave glmWeightedF zinbFit zinbwave
"_PACKAGE"



# FIXME Some of these are redundant with basejump...
#' Parameters
#'
#' @name params
#' @keywords internal
#'
#' @param object Object.
#' @param value Object to assign.
#' @param x Object.
#' @param y Secondary object.
#' @param ... Additional arguments.
#'
#' @param alpha `scalar numeric`. Alpha cutoff (adjusted P value; false
#'   discovery rate).
#' @param aspectRatio `scalar integer`. Aspect ratio.
#' @param bpparam `bpparamClass`. [BiocParallel][] parameter specifying the
#'   back-end to be used for computations. See [BiocParallel::bpparam()] for
#'   details. We recommend one of the following:
#'   - [BiocParallel::bpparam()].
#'   - [BiocParallel::SerialParam()].
#'   - [BiocParallel::MulticoreParam()].
#' [BiocParallel]: https://doi.org/doi:10.18129/B9.bioc.BiocParallel
#' @param color `ggproto`/`ScaleDiscrete` or `NULL`. Desired ggplot2 color
#'   scale. Must supply discrete values. When set to `NULL`, the default ggplot2
#'   color palette will be used. If manual color definitions are desired, we
#'   recommend using [ggplot2::scale_color_manual()].
#'   To set the discrete color palette globally, use
#'   `options(bcbio.discrete.color = scale_color_viridis_d())`.
#' @param dark `boolean`. Plot against a dark background using
#'   [basejump::theme_midnight()].
#' @param dimsUse `integer`. Vector of length 2 that denotes the columns from
#'   the reduced dimension matrix to use for `centerX` and `centerY` column
#'   calculations. Defaults the first and second dimensions.
#' @param dir `string`. Output directory path.
#' @param ensemblRelease `scalar integer`. Ensembl release version (e.g. `90`).
#' @param expression `string`. Calculation to apply. Uses [base::match.arg()]
#'   and defaults to the first argument in the `character` vector.
#' @param fill `ggproto`/`ScaleDiscrete` or `NULL`. Desired ggplot2 fill scale.
#'   Must supply discrete values. When set to `NULL`, the default ggplot2 color
#'   palette will be used. If manual color definitions are desired, we recommend
#'   using [ggplot2::scale_fill_manual()].
#'   To set the discrete fill palette globally, use
#'   `options(bcbio.discrete.fill = scale_fill_viridis_d())`.
#' @param gene2symbol `Gene2Symbol`. Gene-to-symbol mappings.
#' @param genes `character`. Gene identifiers. Must match the rownames of the
#'   object.
#' @param geom `string`. Plot type. Uses [base::match.arg()] and defaults to the
#'   first argument in the `character` vector.
#' @param headerLevel `scalar integer` (`1`-`7`). R Markdown header level.
#' @param interestingGroups `character` or `NULL`. Character vector of
#'   interesting groups. Must be formatted in camel case and intersect with
#'   [sampleData()] colnames.
#' @param label `boolean`. Overlay a cluster identitiy label on the plot.
#' @param labelSize `scalar integer`. Size of the text label.
#' @param legend `boolean`. Include plot legend.
#' @param min `scalar numeric`. Recommended minimum value cutoff.
#' @param max `scalar numeric`. Recommended maximum value cutoff.
#' @param organism `string`. Full Latin organism name (e.g. `"Homo sapiens"`).
#' @param perSample `boolean`. Visualize the distributions per sample?
#' @param pointAlpha `scalar numeric` (`0`-`1`). Alpha transparency level.
#'   Useful when there many cells in the dataset, and some cells can be masked.
#' @param pointsAsNumbers `boolean`. Plot the points as numbers (`TRUE`) or
#'   dots (`FALSE`).
#' @param pointSize `scalar numeric`. Cell point size.
#' @param progress `boolean`. Show progress, using progress bars.
#' @param reducedDim `string`. Name of reduced dimension matrix slotted in
#'   [reducedDims()]. Includes TNSE, UMAP, PCA, for example.
#' @param return `string`. Return type. Uses [base::match.arg()] internally and
#'   defaults to the first argument in the `character` vector.
#' @param title `string` or `NULL`. Plot title.
#' @param trans `string`. Name of the axis scale transformation to apply. See
#'   `help("scale_x_continuous", "ggplot2")` for more information.
#'
#' @return No value.
NULL
