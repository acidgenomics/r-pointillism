#' General Arguments
#'
#' @name general
#' @keywords internal
#'
#' @param object Object.
#' @param value Object to assign.
#' @param x Object.
#' @param y Secondary object.
#' @param ... Additional arguments.
#'
#' @param aspectRatio `scalar integer`. Aspect ratio.
#' @param BPPARAM BiocParallel param. We recommend one of the following:
#'   - [BiocParallel::bpparam()].
#'   - [BiocParallel::SerialParam()].
#'   - [BiocParallel::MulticoreParam()].
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
#' @param expression `string`. Calculation to apply. Uses [match.arg()] and
#'   defaults to the first argument in the `character` vector.
#' @param fill `ggproto`/`ScaleDiscrete` or `NULL`. Desired ggplot2 fill scale.
#'   Must supply discrete values. When set to `NULL`, the default ggplot2 color
#'   palette will be used. If manual color definitions are desired, we recommend
#'   using [ggplot2::scale_fill_manual()].
#'   To set the discrete fill palette globally, use
#'   `options(bcbio.discrete.fill = scale_fill_viridis_d())`.
#' @param genes `character`. Gene identifiers. Must match the rownames of the
#'   object.
#' @param geom `string`. Plot type. Uses [match.arg()] and defaults to the first
#'   argument in the `character` vector.
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
