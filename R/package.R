#' pointillism
#'
#' R package for for single-cell RNA-seq clustering analysis.
#'
#' @keywords internal
"_PACKAGE"



## S4 classes ==================================================================

## Disabled until monocle3 is on Bioconductor.
## > #' @importClassesFrom monocle3 cell_data_set

#' @importClassesFrom AcidGenomes GeneToSymbol
#' @importClassesFrom AcidSingleCell KnownMarkers
#' @importClassesFrom IRanges CompressedSplitDFrameList
#' @importClassesFrom S4Vectors DFrame
#' @importClassesFrom Seurat Seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
NULL



## S4 generics and methods =====================================================

#' @importFrom AcidExperiment sampleNames
#' @importFrom AcidGenerics GeneToSymbol KnownMarkers as.Seurat
#' as.SingleCellExperiment camelCase cellCountsPerCluster cellTypesPerCluster
#' clusters convertGenesToSymbols cpm diffExp diffExpPerCluster findMarkers
#' geometricMean interestingGroups interestingGroups<- leftJoin makeLabel
#' makeNames mapGenesToIds mapGenesToRownames mapGenesToSymbols melt metrics
#' mutateIf plotCellCountsPerCluster plotCellTypesPerCluster plotCounts
#' plotDots plotFeature plotKnownMarkers plotMarker plotPcElbow plotPca
#' plotReducedDim plotStackedBarPlot plotTsne plotTopMarkers plotUmap plotViolin
#' sampleData sampleData<- snakeCase topMarkers uniteInterestingGroups
#' @importFrom BiocGenerics as.data.frame cbind counts counts<- do.call lapply
#' normalize organism organism<- sapply sizeFactors sizeFactors<- sort t table
#' unique unsplit
#' @importFrom BiocParallel bpprogressbar bpprogressbar<-
#' @importFrom Matrix rowMeans rowSums
#' @importFrom S4Vectors as.matrix complete.cases decode head mcols mcols<-
#' metadata metadata<- na.omit split summary tail
#' @importFrom SingleCellExperiment logcounts logcounts<- normcounts
#' normcounts<- reducedDim reducedDim<- reducedDimNames reducedDimNames<-
#' reducedDims
#' @importFrom SummarizedExperiment assay assay<- assayNames assays assays<-
#' colData colData<- rowData rowData<- rowRanges rowRanges<-
#' @importFrom methods coerce show
#' @importFrom pipette import
NULL

#' @importMethodsFrom AcidBase geometricMean
#' @importMethodsFrom AcidExperiment convertGenesToSymbols interestingGroups
#' interestingGroups<- mapGenesToIds mapGenesToRownames mapGenesToSymbols
#' metrics sampleData sampleData<- sampleNames uniteInterestingGroups
#' @importMethodsFrom AcidPlots plotCellCountsPerCluster plotCellTypesPerCluster
#' plotCounts plotDots plotFeature plotKnownMarkers plotMarker plotPca
#' plotReducedDim plotStackedBarPlot plotTsne plotUmap plotViolin
#' @importMethodsFrom AcidPlyr leftJoin melt mutateIf
#' @importMethodsFrom AcidSingleCell cellCountsPerCluster cellTypesPerCluster
#' clusters cpm diffExp diffExpPerCluster findMarkers geometricMean normalize
#' @importMethodsFrom pipette import
#' @importMethodsFrom syntactic camelCase makeLabel makeNames snakeCase
NULL



## S3 generics =================================================================

#' @importFrom Seurat CreateSeuratObject DefaultAssay Idents NormalizeData
#' ScaleData Stdev VariableFeatures
#' @importFrom SeuratObject LayerData Layers
#' @importFrom stats median model.matrix relevel
NULL



## Standard functions ==========================================================

#' @importFrom AcidBase printString showSlotInfo
#' @importFrom AcidCLI abort alert alertInfo alertSuccess alertWarning dl h1 h2
#' separator toInlineString ul verbatim
#' @importFrom AcidExperiment matchInterestingGroups
#' @importFrom AcidGenomes makeGeneToSymbolFromEnsembl
#' @importFrom AcidMarkdown markdownHeader
#' @importFrom AcidPlots .data acid_theme_dark acid_theme_light matchLabels
#' percent pretty_breaks wrap_plots
#' @importFrom BiocParallel bpparam
#' @importFrom IRanges SplitDataFrameList
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom ggplot2 aes element_rect expand_limits facet_grid facet_wrap
#' geom_bar geom_hline geom_jitter geom_line geom_point geom_text geom_violin
#' geom_vline ggplot guide_colorbar guides labs scale_color_gradient
#' scale_color_viridis_c scale_radius scale_x_continuous scale_y_continuous
#' theme theme_set vars
#' @importFrom goalie allAreInLeftOpenRange allAreMatchingRegex
#' allAreNonNegative allArePositive areDisjointSets areIntersectingSets
#' areSetEqual assert bapply false isAlpha isInt hasColnames hasLength
#' hasNames hasRownames hasRows hasValidDimnames hasValidNames isGgscale isAny
#' isCharacter isFlag isHeaderLevel isIntegerish isNonNegative isNumber
#' isPositive isScalar isScalarNumeric isString isSubset quietly
#' requireNamespaces validate
#' @importFrom methods as getMethod is new setAs setClass setMethod
#' setReplaceMethod setValidity slot slot<- validObject .hasSlot
#' @importFrom utils capture.output data packageName packageVersion sessionInfo
NULL
