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
#' @importFrom AcidGenerics as.SingleCellExperiment camelCase
#' @importFrom AcidGenerics cellCountsPerCluster cellTypesPerCluster clusters
#' @importFrom AcidGenerics convertGenesToSymbols cpm diffExp diffExpPerCluster
#' @importFrom AcidGenerics findMarkers geometricMean
#' @importFrom AcidGenerics interestingGroups interestingGroups<-
#' @importFrom AcidGenerics leftJoin makeLabel makeNames
#' @importFrom AcidGenerics mapGenesToIds mapGenesToRownames mapGenesToSymbols
#' @importFrom AcidGenerics melt metrics mutateIf
#' @importFrom AcidGenerics plotCellCountsPerCluster plotCellTypesPerCluster
#' @importFrom AcidGenerics plotCounts plotDots plotFeature plotKnownMarkers
#' @importFrom AcidGenerics plotMarker plotPcElbow plotPca plotReducedDim
#' @importFrom AcidGenerics plotStackedBarPlot plotTsne plotTopMarkers
#' @importFrom AcidGenerics plotUmap plotViolin sampleData sampleData<-
#' @importFrom AcidGenerics snakeCase topMarkers uniteInterestingGroups
#' @importFrom BiocGenerics as.data.frame cbind counts counts<- do.call lapply
#' @importFrom BiocGenerics normalize organism organism<- sapply
#' @importFrom BiocGenerics sizeFactors sizeFactors<- sort t table unique
#' @importFrom BiocGenerics unsplit
#' @importFrom BiocParallel bpprogressbar bpprogressbar<-
#' @importFrom Matrix rowMeans rowSums
#' @importFrom S4Vectors as.matrix complete.cases decode head mcols mcols<-
#' @importFrom S4Vectors metadata metadata<- na.omit split summary tail
#' @importFrom SingleCellExperiment logcounts logcounts<- normcounts
#' @importFrom SingleCellExperiment normcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDim<- reducedDimNames
#' @importFrom SingleCellExperiment reducedDimNames<- reducedDims
#' @importFrom SummarizedExperiment assay assay<- assayNames assays assays<-
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<-
#' @importFrom SummarizedExperiment rowRanges rowRanges<-
#' @importFrom methods coerce show
#' @importFrom pipette import
NULL

#' @importMethodsFrom AcidBase geometricMean
#' @importMethodsFrom AcidExperiment convertGenesToSymbols
#' @importMethodsFrom AcidExperiment interestingGroups interestingGroups<-
#' @importMethodsFrom AcidExperiment mapGenesToIds mapGenesToRownames
#' @importMethodsFrom AcidExperiment mapGenesToSymbols metrics
#' @importMethodsFrom AcidExperiment sampleData sampleData<- sampleNames
#' @importMethodsFrom AcidExperiment uniteInterestingGroups
#' @importMethodsFrom AcidPlots plotCellCountsPerCluster plotCellTypesPerCluster
#' @importMethodsFrom AcidPlots plotCounts plotDots plotFeature plotKnownMarkers
#' @importMethodsFrom AcidPlots plotMarker plotPca plotReducedDim
#' @importMethodsFrom AcidPlots plotStackedBarPlot plotTsne plotUmap plotViolin
#' @importMethodsFrom AcidPlyr leftJoin melt mutateIf
#' @importMethodsFrom AcidSingleCell cellCountsPerCluster cellTypesPerCluster
#' @importMethodsFrom AcidSingleCell clusters cpm diffExp diffExpPerCluster
#' @importMethodsFrom AcidSingleCell findMarkers geometricMean normalize
#' @importMethodsFrom pipette import
#' @importMethodsFrom syntactic camelCase makeLabel makeNames snakeCase
NULL


## S3 generics =================================================================

#' @importFrom Seurat CreateSeuratObject DefaultAssay Idents NormalizeData
#' @importFrom Seurat ScaleData Stdev VariableFeatures
#' @importFrom SeuratObject LayerData Layers
#' @importFrom stats median model.matrix relevel
NULL


## Standard functions ==========================================================

#' @importFrom AcidBase printString showSlotInfo
#' @importFrom AcidCLI abort alert alertInfo alertSuccess alertWarning dl h1 h2
#' @importFrom AcidCLI separator toInlineString ul verbatim
#' @importFrom AcidExperiment matchInterestingGroups
#' @importFrom AcidMarkdown markdownHeader
#' @importFrom AcidPlots .data acid_theme_dark acid_theme_light matchLabels
#' @importFrom AcidPlots percent pretty_breaks wrap_plots
#' @importFrom BiocParallel bpparam
#' @importFrom IRanges SplitDataFrameList
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom ggplot2 aes element_rect expand_limits facet_grid facet_wrap
#' @importFrom ggplot2 geom_bar geom_hline geom_jitter geom_line geom_point
#' @importFrom ggplot2 geom_text geom_violin geom_vline ggplot guide_colorbar
#' @importFrom ggplot2 guides labs scale_color_gradient scale_color_viridis_c
#' @importFrom ggplot2 scale_radius scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 theme theme_set vars
#' @importFrom goalie allAreInLeftOpenRange allAreMatchingRegex
#' @importFrom goalie allAreNonNegative
#' @importFrom goalie allArePositive areDisjointSets areIntersectingSets
#' @importFrom goalie areSetEqual assert bapply false isAlpha isInt hasColnames
#' @importFrom goalie hasLength hasNames hasRownames hasRows hasValidDimnames
#' @importFrom goalie hasValidNames isGgscale isAny isCharacter isFlag
#' @importFrom goalie isHeaderLevel isIntegerish isNonNegative isNumber
#' @importFrom goalie isPositive isScalar isScalarNumeric isString isSubset
#' @importFrom goalie quietly requireNamespaces validate
#' @importFrom methods as getMethod is new setAs setClass setMethod
#' @importFrom methods setReplaceMethod setValidity slot slot<- validObject
#' @importFrom methods .hasSlot
#' @importFrom utils capture.output data packageName packageVersion sessionInfo
NULL
