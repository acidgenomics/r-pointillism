## Conflicts:
## @importFrom Matrix rowMeans rowSums



#' @importClassesFrom Seurat Seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom monocle3 cell_data_set
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom BiocGenerics cbind counts counts<- do.call lapply rowMeans
#'   rowSums sapply
#' @importFrom BiocParallel MulticoreParam SerialParam bpparam bpprogressbar
#'   bpprogressbar<-
#' @importFrom DESeq2 DESeqDataSet DESeq results
#' @importFrom S4Vectors as.data.frame as.matrix complete.cases head mcols
#'   mcols<- metadata metadata<- na.omit split tail
#' @importFrom Seurat CreateSeuratObject Idents as.SingleCellExperiment
#'   as.Seurat GetAssayData Stdev VariableFeatures
#' @importFrom SingleCellExperiment reducedDimNames reducedDims weights
#'   weights<-
#' @importFrom SummarizedExperiment assay assayNames assays assays<- colData
#'   rowData rowRanges rowRanges<-
#' @importFrom acidplots acid_theme_dark acid_theme_light
#' @importFrom basejump camel coerce convertGenesToSymbols decode gene2symbol
#'   import makeNames mapGenesToRownames mapGenesToSymbols markdownHeader
#'   matchArgsToDoCall matchInterestingGroups printString separator showSlotInfo
#'   snake upperCamel
#' @importFrom cowplot plot_grid
#' @importFrom dplyr arrange desc everything filter group_by group_vars
#'   left_join mutate mutate_at mutate_if n pull rename select slice summarize
#'   ungroup vars
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT
#' @importFrom ggplot2 aes element_rect expand_limits facet_grid facet_wrap
#'   geom_bar geom_hline geom_jitter geom_line geom_point geom_text geom_violin
#'   geom_vline ggplot guide_colorbar guides labs scale_color_gradient
#'   scale_color_viridis_c scale_radius scale_x_continuous scale_y_continuous
#'   theme theme_set
#' @importFrom goalie allAreInLeftOpenRange allAreNonNegative allArePositive
#'   areDisjointSets areIntersectingSets assert false isAlpha isInt hasLength
#'   hasNames hasRownames hasRows isGGScale isAny isCharacter isFlag
#'   isHeaderLevel isIntegerish isNonEmpty isNumber isPositive isScalar isString
#'   isSubset validate
#' @importFrom magrittr %>% set_names set_rownames
#' @importFrom methods as getMethod is new setAs setClass show slot slot<-
#'   validObject
#' @importFrom pbapply pblapply
#' @importFrom rlang !! !!! := sym syms
#' @importFrom scales percent pretty_breaks
#' @importFrom sessioninfo session_info
#' @importFrom stats median model.matrix relevel
#' @importFrom tibble as_tibble column_to_rownames remove_rownames tibble
#' @importFrom tidyr gather
#' @importFrom utils capture.output globalVariables packageVersion
NULL
