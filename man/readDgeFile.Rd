% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fileParsing.R
\name{readDgeFile}
\alias{readDgeFile}
\title{Read in DGE matrix}
\usage{
readDgeFile(dgeMatrixFile, cell_features = NULL)
}
\arguments{
\item{dgeMatrixFile}{The file to parse.  If a directory, then the expected
10x MTX files are read.}

\item{cell_features}{A cell features data frame with the cell barcodes to
match the DGE matrix.  If supplied, the DGE matrix will be reordered and
filtered to match the cell_features barcodes.}
}
\value{
A matrix of the DGE matrix with (optionally) the cell barcodes in the
same order as the cell_features
}
\description{
Read in a dense DGE matrix in McCarroll lab format (cell barcodes in columns,
genes in rows) or read the expected 10x MTX files from a directory.
}
\details{
Validates and reorders DGE matrix to match the cell_features order.
}
\examples{
# Create a temporary directory for test data
temp_dir <- tempdir()

# Write example dense DGE matrix to the temporary directory
dge_file <- writeExampleDenseDGE(temp_dir)

# Write example cell features to the temporary directory
cell_features_file=writeExampleCellFeatures(temp_dir)

# Load example cell features
cell_features <- read.table(cell_features_file, sep = "\t", header = TRUE,
stringsAsFactors = FALSE)

# Read the DGE matrix and reorder it to match the cell_features barcodes
dge_matrix <- readDgeFile(dge_file, cell_features)

# Display output summary
dim(dge_matrix)  # Check dimensions of the parsed DGE matrix
head(dge_matrix[, 1:5])  # Preview the first few columns

# Read the DGE matrix and without any cell feature dependency
dge_matrix <- readDgeFile(dge_file, cell_features=NULL)

# Display output summary
dim(dge_matrix)  # Check dimensions of the parsed DGE matrix
# Cell barcodes have a different sort order.
head(dge_matrix[, 1:5])  # Preview the first few columns
}
