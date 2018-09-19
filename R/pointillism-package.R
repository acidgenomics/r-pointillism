#' pointillism
#'
#' R package for for single-cell RNA-seq clustering analysis.
#'
#' @name pointillism-package
#' @docType package
#'
#' @importClassesFrom Seurat seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @importFrom assertive.base assert_are_identical assert_all_are_not_na
#' @importFrom assertive.numbers assert_all_are_in_left_open_range
#'   assert_all_are_positive
#' @importFrom assertive.properties assert_has_no_duplicates assert_has_names
#'   assert_has_rows assert_is_non_empty assert_is_of_length
#' @importFrom assertive.sets assert_are_disjoint_sets
#'   assert_are_intersecting_sets assert_is_subset
#' @importFrom assertive.types assert_is_a_bool assert_is_a_number
#'   assert_is_a_string assert_is_all_of assert_is_any_of assert_is_character
#'   assert_is_data.frame assert_is_factor assert_is_list assert_is_matrix
#'   assert_is_numeric assert_is_tbl_df is_a_string
#' @importFrom basejump assertHasRownames assertIsAHeaderLevel
#'   assertIsAStringOrNULL assertIsAnImplicitInteger
#'   assertIsColorScaleContinuousOrNULL assertIsColorScaleDiscreteOrNULL
#'   assertIsFillScaleDiscreteOrNULL assertIsImplicitInteger camel coerce
#'   convertGenesToSymbols gene2symbol import mapGenesToRownames
#'   mapGenesToSymbols markdownHeader matchInterestingGroups printString
#'   theme_midnight uniteInterestingGroups upperCamel
#' @importFrom bcbioSingleCell filterCells
#' @importFrom BiocGenerics cbind counts counts<- do.call
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
#' @importFrom magrittr %>% set_names set_rownames
#' @importFrom Matrix rowMeans rowSums
#' @importFrom methods as getMethod is setAs show slot slot<- validObject
#' @importFrom pbapply pblapply
#' @importFrom rlang !! !!! sym syms
#' @importFrom S4Vectors as.data.frame as.matrix
#'   complete.cases mcols mcols<- metadata metadata<- na.omit
#' @importFrom scales percent pretty_breaks
#' @importFrom Seurat as.SingleCellExperiment CreateSeuratObject
#' @importFrom stats median model.matrix relevel
#' @importFrom SummarizedExperiment assay assayNames assays assays<- colData
#'   rowData rowRanges rowRanges<-
#' @importFrom tibble remove_rownames tibble
#' @importFrom tidyr gather
#' @importFrom utils capture.output globalVariables packageVersion
#' @importFrom zinbwave glmWeightedF zinbFit zinbwave
NULL
