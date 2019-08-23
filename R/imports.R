## Conflicts:
## #' @importFrom Matrix rowMeans rowSums



#' @importClassesFrom Seurat Seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom monocle3 cell_data_set
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom BiocGenerics cbind counts counts<- do.call lapply rowMeans
#'   rowSums sapply sizeFactors sizeFactors<-
#' @importFrom BiocParallel MulticoreParam SerialParam bpparam bpprogressbar
#'   bpprogressbar<-
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colSums2 rowSums2
#' @importFrom DESeq2 DESeqDataSet DESeq estimateSizeFactorsForMatrix results
#' @importFrom S4Vectors DataFrame as.data.frame as.matrix complete.cases head
#'   mcols mcols<- metadata metadata<- na.omit split tail
#' @importFrom Seurat CreateSeuratObject Idents as.SingleCellExperiment
#'   as.Seurat DefaultAssay GetAssayData NormalizeData Stdev VariableFeatures
#' @importFrom SingleCellExperiment logcounts logcounts<- normcounts
#'   normcounts<- reducedDim reducedDim<- reducedDimNames reducedDims
#'   sizeFactors sizeFactors<- spikeNames
#' @importFrom SummarizedExperiment assay assay<- assayNames assays assays<-
#'   colData colData<- rowData rowData<- rowRanges rowRanges<-
#' @importFrom acidplots acid_theme_dark acid_theme_light
#' @importFrom basejump Gene2Symbol camelCase coerce convertGenesToSymbols
#'   decode estimateSizeFactors gene2symbol geometricMean import
#'   interestingGroups interestingGroups<- makeNames mapGenesToIDs
#'   mapGenesToRownames mapGenesToSymbols markdownHeader matchArgsToDoCall
#'   matchInterestingGroups metrics organism organism<- printString sampleData
#'   sampleData<- sampleNames separator showSlotInfo snakeCase upperCamel
#' @importFrom cowplot plot_grid
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT
#' @importFrom ggplot2 aes element_rect expand_limits facet_grid facet_wrap
#'   geom_bar geom_hline geom_jitter geom_line geom_point geom_text geom_violin
#'   geom_vline ggplot guide_colorbar guides labs scale_color_gradient
#'   scale_color_viridis_c scale_radius scale_x_continuous scale_y_continuous
#'   theme theme_set
#' @importFrom goalie allAreInLeftOpenRange allAreNonNegative allArePositive
#'   areDisjointSets areIntersectingSets assert bapply false isAlpha isInt
#'   hasColnames hasLength hasNames hasRownames hasRows hasValidDimnames
#'   hasValidNames isGGScale isAny isCharacter isFlag isHeaderLevel isIntegerish
#'   isNonEmpty isNonNegative isNumber isPositive isScalar isScalarNumeric
#'   isString isSubset validate
#' @importFrom methods as getMethod is new setAs setClass setMethod
#'   setReplaceMethod setValidity show slot slot<- validObject .hasSlot
#' @importFrom rlang !! !!! := sym syms
#' @importFrom scales percent pretty_breaks
#' @importFrom scater calculateCPM normalizeSCE
#' @importFrom sessioninfo session_info
#' @importFrom stats median model.matrix relevel
#' @importFrom utils capture.output globalVariables packageVersion
#'
#'
#'
#' @importFrom dplyr arrange desc everything filter group_by group_vars
#'   left_join mutate mutate_at mutate_if n pull rename select slice summarize
#'   ungroup vars
#' @importFrom magrittr %>% set_names set_rownames
#' @importFrom tidyr gather
#' @importFrom tibble as_tibble column_to_rownames remove_rownames
#'   rownames_to_column tibble
NULL
