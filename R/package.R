#' pointillism
#'
#' R package for for single-cell RNA-seq clustering analysis.
#'
#' @keywords internal
#'
#' @importClassesFrom basejump SingleCellExperiment
#' @importClassesFrom Seurat Seurat
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom BiocParallel bpparam bpprogressbar bpprogressbar<-
#' @importFrom DESeq2 DESeqDataSet DESeq estimateSizeFactorsForMatrix results
#' @importFrom Seurat CreateSeuratObject DefaultAssay GetAssayData Idents
#'   NormalizeData Stdev VariableFeatures as.SingleCellExperiment as.Seurat
#' @importFrom AcidPlots !! !!! acid_theme_dark acid_theme_light matchLabels
#'   percent plot_grid pretty_breaks sym syms
#' @importFrom basejump DataFrame Gene2Symbol SingleCellExperiment
#'   SplitDataFrameList alert alertInfo alertSuccess alertWarning as.data.frame
#'   as.matrix assay assay<- assayNames assays assays<- camelCase capture.output
#'   cbind coerce colData colData<- complete.cases convertGenesToSymbols counts
#'   counts<- data decode do.call estimateSizeFactors geometricMean h1 h2 head
#'   import interestingGroups interestingGroups<- lapply leftJoin logcounts
#'   logcounts<- makeGene2SymbolFromEnsembl makeNames mapGenesToIDs
#'   mapGenesToRownames mapGenesToSymbols markdownHeader matchInterestingGroups
#'   mcols mcols<- median melt metadata metadata<- metrics model.matrix mutateIf
#'   normcounts normcounts<- organism organism<- na.omit packageName
#'   packageVersion printString reducedDim reducedDim<- reducedDimNames
#'   reducedDims relevel rowData rowData<- rowMeans rowRanges rowRanges<-
#'   rowSums sampleData sampleData<- sampleNames sapply separator session_info
#'   showSlotInfo sizeFactors sizeFactors<- snakeCase split t table tail ul
#'   uniteInterestingGroups unsplit verbatim
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
