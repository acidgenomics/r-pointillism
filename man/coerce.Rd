% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coerce-methods.R
\name{coerce}
\alias{coerce}
\alias{coerce,Seurat,SingleCellExperiment-method}
\alias{coerce,Seurat,RangedSummarizedExperiment-method}
\alias{coerce,Seurat,SummarizedExperiment-method}
\alias{coerce,SingleCellExperiment,Seurat-method}
\title{Force an object to belong to a class}
\description{
Force an object to belong to a class
}
\note{
Updated 2021-03-03.
}
\section{\code{Seurat} to \code{SingleCellExperiment}}{

S4 coercion support for creating a \code{SingleCellExperiment} from a \code{Seurat}
class object. The \href{https://satijalab.org/Seurat/faq}{Seurat FAQ page}
explains the \code{Seurat} S4 class structure in detail. Internally, this method
improves the basic \code{Seurat::as.SingleCellExperiment} S3 coercion method,
including the \code{object@scale.data} matrix, and will keep track of stashed
\code{rowRanges} and \code{metadata} if the \code{Seurat} object was originally created
from a \code{SingleCellExperiment} (i.e. from the bcbioSingleCell package).
}

\section{\code{Seurat} to \code{RangedSummarizedExperiment}}{

S4 coercion support for creating a \code{RangedSummarizedExperiment} from a
\code{Seurat} class object.
}

\section{\code{Seurat} to \code{SummarizedExperiment}}{

S4 coercion support for creating a \code{SummarizedExperiment} from a \code{Seurat}
class object.
}

\section{\code{SingleCellExperiment} to \code{Seurat}}{

Interally \code{Seurat::CreateSeuratObject} is called without applying any
additional filtering cutoffs, since we have already defined them during our
quality control analysis. Here we are passing the raw gene-level counts of
the filtered cells into a new \code{Seurat} class object. Use
\code{convertGenesToSymbols} to convert gene IDs to names (symbols).
}

\examples{
data(Seurat, SingleCellExperiment, package = "AcidTest")

## SingleCellExperiment to Seurat ====
x <- as(SingleCellExperiment, "Seurat")
class(x)
print(x)

## Seurat to SingleCellExperiment ====
x <- as(Seurat, "SingleCellExperiment")
print(x)
}
\seealso{
\itemize{
\item \code{Seurat::CreateSeuratObject()}.
\item \code{Seurat::as.Seurat()}.
\item \code{Seurat::as.SingleCellExperiment()}.
}
}
