#' pointillism
#'
#' R package for for single-cell RNA-seq clustering analysis.
#'
#' @keywords internal
#'
#' @importClassesFrom Seurat Seurat
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom DESeq2 DESeqDataSet DESeq estimateSizeFactorsForMatrix results
#' @importFrom S4Vectors DataFrame as.data.frame as.matrix cbind complete.cases
#'   do.call head lapply mcols mcols<- metadata metadata<- na.omit sapply split
#'   tail
#' @importFrom Seurat CreateSeuratObject DefaultAssay GetAssayData Idents
#'   NormalizeData Stdev VariableFeatures as.SingleCellExperiment as.Seurat
#' @importFrom SingleCellExperiment logcounts logcounts<- normcounts
#'   normcounts<- reducedDim reducedDim<- reducedDimNames reducedDims
#'   sizeFactors sizeFactors<-
#' @importFrom AcidPlots !! !!! acid_theme_dark acid_theme_light matchLabels
#'   plot_grid sym syms
#' @importFrom basejump Gene2Symbol SingleCellExperiment SplitDataFrameList
#'   assay assay<- assayNames assays assays<- camelCase capture.output coerce
#'   colData colData<- convertGenesToSymbols counts counts<- decode
#'   estimateSizeFactors geometricMean import interestingGroups
#'   interestingGroups<- leftJoin makeGene2SymbolFromEnsembl makeNames
#'   mapGenesToIDs mapGenesToRownames mapGenesToSymbols markdownHeader
#'   matchInterestingGroups melt metrics mutateIf organism organism<-
#'   packageName packageVersion printString rowData rowData<- rowMeans rowRanges
#'   rowRanges<- rowSums sampleData sampleData<- sampleNames separator
#'   session_info showSlotInfo sizeFactors sizeFactors<- snakeCase table
#'   uniteInterestingGroups unsplit
#' @importFrom cli cli_alert cli_alert_info cli_alert_success cli_alert_warning
#'   cli_div cli_dl cli_end cli_h1 cli_h2 cli_text cli_ul
#' @importFrom dplyr group_by n summarize
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT
#' @importFrom ggplot2 aes element_rect expand_limits facet_grid facet_wrap
#'   geom_bar geom_hline geom_jitter geom_line geom_point geom_text geom_violin
#'   geom_vline ggplot guide_colorbar guides labs scale_color_gradient
#'   scale_color_viridis_c scale_radius scale_x_continuous scale_y_continuous
#'   theme theme_set vars
#' @importFrom goalie allAreInLeftOpenRange allAreMatchingRegex
#'   allAreNonNegative allArePositive areDisjointSets areIntersectingSets
#'   areSetEqual assert bapply false isAlpha isInt hasColnames hasLength
#'   hasNames hasRownames hasRows hasValidDimnames hasValidNames isGGScale isAny
#'   isCharacter isFlag isHeaderLevel isIntegerish isNonNegative isNumber
#'   isPositive isScalar isScalarNumeric isString isSubset validate
#' @importFrom methods as getMethod is new setAs setClass setMethod
#'   setReplaceMethod setValidity show slot slot<- validObject .hasSlot
#' @importFrom scater calculateCPM normalizeCounts
"_PACKAGE"



## Disabled until monocle3 is on Bioconductor.
## > #' @importClassesFrom monocle3 cell_data_set

## Conflicts (FIXME THIS IS OK):
## > #' @importFrom Matrix rowMeans rowSums



## FIXME Take these out?
#' @importFrom BiocParallel MulticoreParam SerialParam bpparam bpprogressbar
#'   bpprogressbar<-
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colSums2 rowSums2
#' @importFrom scales percent pretty_breaks
#' @importFrom stats median model.matrix relevel
#' @importFrom utils data
NULL
