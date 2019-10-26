## Conflicts:
## #' @importFrom Matrix rowMeans rowSums



#' @importClassesFrom Seurat Seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom monocle3 cell_data_set
NULL



#' @importMethodsFrom basejump coerce
NULL



#' @importFrom BiocGenerics counts counts<- rowMeans rowSums sizeFactors
#'   sizeFactors<- table
#' @importFrom BiocParallel MulticoreParam SerialParam bpparam bpprogressbar
#'   bpprogressbar<-
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colSums2 rowSums2
#' @importFrom DESeq2 DESeqDataSet DESeq estimateSizeFactorsForMatrix results
#' @importFrom IRanges SplitDataFrameList unsplit
#' @importFrom S4Vectors DataFrame as.data.frame as.matrix cbind complete.cases
#'   do.call head lapply mcols mcols<- metadata metadata<- na.omit sapply split
#'   tail
#' @importFrom Seurat CreateSeuratObject DefaultAssay GetAssayData Idents
#'   NormalizeData Stdev VariableFeatures as.SingleCellExperiment as.Seurat
#' @importFrom SingleCellExperiment logcounts logcounts<- normcounts
#'   normcounts<- reducedDim reducedDim<- reducedDimNames reducedDims
#'   sizeFactors sizeFactors<- spikeNames
#' @importFrom SummarizedExperiment assay assay<- assayNames assays assays<-
#'   colData colData<- rowData rowData<- rowRanges rowRanges<-
#' @importFrom acidplots acid_theme_dark acid_theme_light
#' @importFrom basejump Gene2Symbol camelCase coerce convertGenesToSymbols
#'   decode estimateSizeFactors geometricMean import interestingGroups
#'   interestingGroups<- leftJoin makeNames mapGenesToIDs mapGenesToRownames
#'   mapGenesToSymbols markdownHeader matchArgsToDoCall matchInterestingGroups
#'   melt metrics mutateIf organism organism<- printString sampleData
#'   sampleData<- sampleNames separator showSlotInfo snakeCase
#'   uniteInterestingGroups
#' @importFrom cowplot plot_grid
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
#' @importFrom rlang !! sym
#' @importFrom scales percent pretty_breaks
#' @importFrom scater calculateCPM normalizeSCE
#' @importFrom sessioninfo session_info
#' @importFrom stats median model.matrix relevel
#' @importFrom utils capture.output globalVariables packageVersion
NULL
