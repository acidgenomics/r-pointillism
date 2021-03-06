% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/KnownMarkers-methods.R
\name{KnownMarkers}
\alias{KnownMarkers}
\alias{KnownMarkers,SeuratMarkersPerCluster,CellTypeMarkers-method}
\title{Known markers}
\usage{
KnownMarkers(markers, known, ...)

\S4method{KnownMarkers}{SeuratMarkersPerCluster,CellTypeMarkers}(markers, known, promiscuousThreshold = 0L)
}
\arguments{
\item{markers}{\code{SeuratMarkers} or \code{SeuratMarkersPerCluster}.}

\item{known}{\code{CellTypeMarkers}.
Grouped by \code{cellType} column. Known markers \code{data.frame} imported by
\code{readCellTypeMarkers} or pulled from internal cell cycle markers data.}

\item{...}{Additional arguments.}

\item{promiscuousThreshold}{\code{integer(1)}.
Minimum number of clusters required to consider a gene marker promiscuous.
Set to \code{0} to disable promiscuous marker filtering.}
}
\value{
\code{KnownMarkers}.
}
\description{
Known markers
}
\note{
Both the \code{markers} and \code{known} objects must contain Ensembl gene
identifiers in the \code{geneId} column. We must avoid any matching operations
based on the gene names, since these change often and can mismatch
easily.

Updated 2021-03-03.
}
\examples{
data(cellTypeMarkersList, seuratAllMarkers)

## SeuratMarkersPerCluster ====
markers <- seuratAllMarkers
known <- cellTypeMarkersList[["homoSapiens"]]
x <- KnownMarkers(markers = markers, known = known)
summary(x)
}
