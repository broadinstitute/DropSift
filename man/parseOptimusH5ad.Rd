% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parseOptimus.R
\name{parseOptimusH5ad}
\alias{parseOptimusH5ad}
\title{Parse an Optimus pipeline generated h5ad file.}
\usage{
parseOptimusH5ad(h5ad_file, min_transcripts = 20)
}
\arguments{
\item{h5ad_file}{Path to the h5ad file}

\item{min_transcripts}{Minimum number of transcripts required for a cell barcode to be included in the output.}
}
\value{
A list(dge, cell_features) containing the sparse matrix with genes in rows and cell barcodes in columns,
and a dataframe of cell features needed for nuclei selection.
}
\description{
This function parses the input, extracts the expression data as a sparse matrix with
cells in columns and gens in rows.  The cell level metadata needed for nuclei selection
is calculated / extracted from the obs dataframe.
}
