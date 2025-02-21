# if (!require('BiocManager', quietly = TRUE))
# install.packages('BiocManager') BiocManager::install('rhdf5')


# library(rhdf5) library(Matrix)


#' Parse an H5AD file and extract the expression matrix and cell features (obs
#' dataframe)
#'
#' This function tests if the incoming data is organized with genes in rows and
#' cells in columns or the transpose. Output data is always genes in rows and
#' cells in columns (Seurat/DropSeq format).
#'
#' If the file is gzipped (.h5ad.gz), it is automatically decompressed to a
#' temporary file before reading.
#'
#' @param h5ad_file Path to the h5ad file (can be .h5ad or .h5ad.gz)
#' @param expression_matrix_path Path to the raw counts expression matrix in the
#'   h5ad file
#' @param gene_id_path Path to the gene IDs in the h5ad file.
#' @param cell_id_path Path to the cell IDs in the h5ad file.
#' @return A list containing the sparse matrix with genes in rows and cell
#'   barcodes in columns, and a dataframe of cell features if the obs path
#'   exists in the h5ad file.
#' @import rhdf5 Matrix R.utils
#' @export
#' @examples
#' # Load the example h5ad file
#' h5ad_file <- system.file("extdata", "adata_example.h5ad.gz", package =
#' "DropSift")
#' r <- parseH5ad(h5ad_file)
#' dge <- r$dge
#' cell_features <- r$cell_features
#' # Display output summary
#'
#' # Check dimensions of the parsed expression matrix
#' dim(dge)
#' r$dge # The example DGE is 5 cells and 5 genes.
#'
#' # This cell feature set is the raw outputs from Optimus which are not yet
#' # read to be used by DropSift.  See parseOptimusH5ad.
#' head(cell_features)
parseH5ad <- function(h5ad_file, expression_matrix_path = "X",
    gene_id_path = "/var/ensembl_ids",
    cell_id_path = "/obs/CellID") {

    if (!file.exists(h5ad_file))
        stop("The file [", h5ad_file, "] does not exist.")

    # Detect if the file is gzipped
    is_gzipped <- grepl("\\.gz$", h5ad_file)
    # If gzipped, decompress to a temporary file
    if (is_gzipped) {
        temp_file <- tempfile(fileext = ".h5ad")
        R.utils::gunzip(h5ad_file, destname = temp_file, remove = FALSE)
        h5ad_file <- temp_file  # Use decompressed file
    }

    # Load metadata and expression matrix
    metadata <- load_h5ad_metadata(h5ad_file, gene_id_path, cell_id_path)
    expr_data <- load_h5ad_expression(h5ad_file, expression_matrix_path)

    # Determine matrix orientation
    orientation <- get_matrix_orientation(
        expr_data$csr_indices, expr_data$csr_indptr,
        length(metadata$ensembl_ids), length(metadata$cell_names)
    )
    # Create the sparse matrix
    sparse_matrix <- sparseMatrix(
        i = orientation$row_indices,
        j = orientation$col_indices,
        x = expr_data$csr_data,
        dims = c(orientation$matrix_rows, orientation$matrix_cols)
    )
    # Transpose if needed
    if (orientation$transpose_needed) {
        sparse_matrix <- t(sparse_matrix)
    }
    # Assign row and column names
    rownames(sparse_matrix) <- metadata$ensembl_ids
    colnames(sparse_matrix) <- metadata$cell_names
    # Load cell features
    obs_df <- load_h5ad_obs(h5ad_file)
    # Clean up temporary file if decompressed
    if (is_gzipped) {
        unlink(temp_file)
    }

    list(dge = sparse_matrix, cell_features = obs_df)
}


#' Load gene and cell names from an h5ad file
#' @noRd
load_h5ad_metadata <- function(h5ad_file, gene_id_path, cell_id_path) {
    ensembl_ids <- as.character(h5read(h5ad_file, gene_id_path))
    cell_names <- as.character(h5read(h5ad_file, cell_id_path))
    list(ensembl_ids = ensembl_ids, cell_names = cell_names)
}

#' Load expression matrix and convert to CSR format
#' @noRd
load_h5ad_expression <- function(h5ad_file, expression_matrix_path) {
    expression_matrix <- rhdf5::h5read(h5ad_file, expression_matrix_path)
    list(
        csr_data = as.numeric(expression_matrix$data),
        csr_indices = expression_matrix$indices,
        csr_indptr = expression_matrix$indptr
    )
}

#' Determine matrix dimensions and check orientation
#' @noRd
get_matrix_orientation <- function(csr_indices,
    csr_indptr, num_genes, num_cells) {

    matrix_rows <- length(csr_indptr) - 1
    matrix_cols <- max(csr_indices) + 1
    transpose_needed <- FALSE

    if (matrix_rows == num_genes && matrix_cols == num_cells) {
        row_indices <- rep(seq_len(matrix_rows), diff(csr_indptr))
        col_indices <- csr_indices + 1  # Adjust for 0-based indexing
    } else if (matrix_rows == num_cells && matrix_cols == num_genes) {
        row_indices <- rep(seq_len(matrix_rows), diff(csr_indptr))
        col_indices <- csr_indices + 1
        transpose_needed <- TRUE
    } else {
        stop("Matrix dimensions do not match metadata dimensions.")
    }

    list(
        matrix_rows = matrix_rows,
        matrix_cols = matrix_cols,
        row_indices = row_indices,
        col_indices = col_indices,
        transpose_needed = transpose_needed
    )
}

#' Load optional `obs` data
#' @noRd
load_h5ad_obs <- function(h5ad_file) {
    obs_path_exists <- any(h5ls(h5ad_file)$name == "obs")
    if (obs_path_exists) {
        as.data.frame(h5read(h5ad_file, "/obs"))
    } else {
        NULL
    }
}





#' Parse an Optimus pipeline generated h5ad file.
#'
#' This function parses the input, extracts the expression data as a sparse
#' matrix with cells in columns and gens in rows.  The cell level metadata
#' needed for nuclei selection is calculated / extracted from the obs dataframe.
#'
#' @param h5ad_file Path to the h5ad file
#' @param min_transcripts Minimum number of transcripts required for a cell
#'   barcode to be included in the output.
#' @return A list(dge, cell_features) containing the sparse matrix with genes in
#'   rows and cell barcodes in columns, and a dataframe of cell features needed
#'   for nuclei selection.
#' @export
#' @examples
#' # Load the example h5ad file
#' h5ad_file <- system.file("extdata", "adata_example.h5ad.gz", package =
#' "DropSift")
#' r <- parseOptimusH5ad(h5ad_file, min_transcripts=0)
#' dge <- r$dge
#' cell_features <- r$cell_features
#' # Display output summary
#' dim(dge)  # Check dimensions of the parsed expression matrix
#' r$dge # The example DGE is 5 cells and 5 genes
#' #This cell feature set is processed to include the features needed for
#' DropSift.
#' head(cell_features)  # Preview the cell features dataframe.
parseOptimusH5ad <- function(h5ad_file, min_transcripts = 20) {
    log_info(paste("Reading expression from Optimus h5ad file [", h5ad_file,
        "]", sep = ""))
    r <- parseH5ad(h5ad_file)
    r <- prepareOptimusDataForNucleiSelection(sparse_matrix = r$dge,
        obs_df = r$cell_features, min_transcripts = min_transcripts)
    return(r)
}


# updates the obs dataframe with the % intronic, % mt, and number of
# transcripts and cell_barcodes columns.  formats to the expectations
# of the SVM expected column names - cell_barcode, num_transcripts,
# num_reads, pct_intronic, pct_mt
prepareOptimusDataForNucleiSelection <- function(sparse_matrix, obs_df,
    min_transcripts = 20) {
    # calculate the % intronic
    obs_df$pct_intronic <- (obs_df$reads_mapped_intronic +
        obs_df$reads_mapped_intronic_as)/obs_df$reads_mapped_uniquely

    # calculate the % mt
    obs_df$pct_mt <-
        (obs_df$reads_mapped_mitochondrial/obs_df$reads_mapped_uniquely)

    # calculate the number of transcripts
    obs_df$num_transcripts <- colSums(sparse_matrix)

    df <- data.frame(cell_barcode = obs_df$CellID,
        num_transcripts = obs_df$num_transcripts,
        num_reads = obs_df$reads_mapped_uniquely,
        pct_intronic = obs_df$pct_intronic,
        pct_mt = obs_df$pct_mt)

    if ("star_IsCell" %in% colnames(obs_df)) {
        df$IsCell <- obs_df$star_IsCell
    }

    # filter to at least some minimum number of transcripts.
    idx <- which(df$num_transcripts >= min_transcripts)
    df <- df[idx, ]

    # filter sparse matrix to same barcodes in same order.
    dge <- sparse_matrix[, df$cell_barcode]

    result <- list(dge = dge, cell_features = df)
    return(result)
}

getOptimusGeneSymbols <- function(h5ad_file) {
    # Get the gene names.  Gene names can be repeated.  Read gene
    # indices and corresponding gene names
    gene_indices <- h5read(h5ad_file, "/var/Gene")  # Gene indices
    gene_names <- h5read(h5ad_file, "/var/__categories/Gene")

    # Map indices to names (if necessary) If the indices are numeric,
    # they map directly to the gene names
    gene_names_full <- gene_names[gene_indices + 1]
    return(as.character(gene_names_full))
}
