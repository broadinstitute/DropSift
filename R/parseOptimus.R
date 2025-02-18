# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("rhdf5")


#library(rhdf5)
#library(Matrix)

#' Parse an 5had file and extract the expression matrix and cell features (obs dataframe)
#'
#' This function tests if the incoming data is organized with genes in rows and cells in columns
#' or the transpose.  Output data is always genes in rows and cells in columns (Seurat/DropSeq format)
#' @param h5ad_file Path to the h5ad file
#' @param expression_matrix_path Path to the raw counts expression matrix in the h5ad file
#' @param gene_id_path Path to the gene IDs in the h5ad file.
#' @param cell_id_path Path to the cell IDs in the h5ad file.
#' @return A list containing the sparse matrix with genes in rows and cell barcodes in columns, and a dataframe of cell features if
#' the /obs path exists in the h5ad file.
#' @import rhdf5 Matrix
#' @export
parseH5ad <- function(h5ad_file, expression_matrix_path = "X",
                      gene_id_path = "/var/ensembl_ids", cell_id_path = "/obs/CellID") {

    if (!file.exists(h5ad_file)) {
        stop("The file [", h5ad_file, "] does not exist.")
    }

    # Load the cell and gene names
    # these need to be characters and not arrays for later on when
    # tranformed into a seurat object.
    ensembl_ids <- as.character(h5read(h5ad_file, gene_id_path))
    cell_names <- as.character(h5read(h5ad_file, cell_id_path))

    # Load the expression matrix (CSR format)
    expression_matrix <- h5read(h5ad_file, expression_matrix_path)
    csr_data <- as.numeric(expression_matrix$data)  # Non-zero values
    csr_indices <- expression_matrix$indices        # Column indices
    csr_indptr <- expression_matrix$indptr          # Row pointer

    # Determine dimensions
    num_genes <- length(ensembl_ids)
    num_cells <- length(cell_names)

    # Check matrix orientation
    matrix_rows <- length(csr_indptr) - 1
    matrix_cols <- max(csr_indices) + 1
    transpose_needed <- FALSE

    if (matrix_rows == num_genes && matrix_cols == num_cells) {
        # Orientation matches: genes in rows, cells in columns
        row_indices <- rep(seq_len(matrix_rows), diff(csr_indptr))
        col_indices <- csr_indices + 1  # Adjust for 0-based indexing
    } else if (matrix_rows == num_cells && matrix_cols == num_genes) {
        # Orientation is transposed: cells in rows, genes in columns
        row_indices <- rep(seq_len(matrix_rows), diff(csr_indptr))
        col_indices <- csr_indices + 1  # Adjust for 0-based indexing
        transpose_needed <- TRUE
    } else {
        stop("Matrix dimensions do not match metadata dimensions.")
    }

    # Create the sparse matrix
    sparse_matrix <- sparseMatrix(
        i = row_indices,
        j = col_indices,
        x = csr_data,
        dims = c(matrix_rows, matrix_cols)
    )

    # Transpose the matrix if needed to ensure genes in rows and cells in columns
    if (transpose_needed) {
        sparse_matrix <- t(sparse_matrix)
    }

    # Ensure row and column names are consistent with orientation
    rownames(sparse_matrix) <- ensembl_ids
    colnames(sparse_matrix) <- cell_names

    # Check if `/obs` exists
    obs_path_exists <- any(h5ls(h5ad_file)$name == "obs")

    # Load `obs` data if it exists
    if (obs_path_exists) {
        obs_data <- h5read(h5ad_file, "/obs")
        obs_df <- as.data.frame(obs_data)
    } else {
        obs_df <- NULL
    }

    result <- list(dge = sparse_matrix, cell_features = obs_df)

    return(result)
}

#' Parse an Optimus pipeline generated h5ad file.
#'
#' This function parses the input, extracts the expression data as a sparse matrix with
#' cells in columns and gens in rows.  The cell level metadata needed for nuclei selection
#' is calculated / extracted from the obs dataframe.
#'
#' @param h5ad_file Path to the h5ad file
#' @param min_transcripts Minimum number of transcripts required for a cell barcode to be included in the output.
#' @return A list(dge, cell_features) containing the sparse matrix with genes in rows and cell barcodes in columns,
#' and a dataframe of cell features needed for nuclei selection.
#' @export
parseOptimusH5ad<-function (h5ad_file, min_transcripts=20) {
    log_info(paste("Reading expression from Optimus h5ad file [", h5ad_file, "]", sep=""))
    r <- parseH5ad(h5ad_file)
    r <- prepareOptimusDataForNucleiSelection(sparse_matrix=r$dge, obs_df=r$cell_features, min_transcripts=min_transcripts)
    return (r)
}


#updates the obs dataframe with the % intronic, % mt, and number of transcripts and cell_barcodes columns.
#formats to the expectations of the SVM expected column names - cell_barcode, num_transcripts, num_reads, pct_intronic, pct_mt
prepareOptimusDataForNucleiSelection<-function(sparse_matrix, obs_df, min_transcripts=20) {
    #calculate the % intronic
    obs_df$pct_intronic <- (obs_df$reads_mapped_intronic+obs_df$reads_mapped_intronic_as)/obs_df$reads_mapped_uniquely
    #calculate the % mt
    obs_df$pct_mt <- (obs_df$reads_mapped_mitochondrial/obs_df$reads_mapped_uniquely)
    #calculate the number of transcripts
    obs_df$num_transcripts <- colSums(sparse_matrix)

    df<- data.frame(cell_barcode=obs_df$CellID,
                    num_transcripts=obs_df$num_transcripts,
                    num_reads=obs_df$reads_mapped_uniquely,
                    pct_intronic=obs_df$pct_intronic,
                    pct_mt=obs_df$pct_mt)

    if ("star_IsCell" %in% colnames(obs_df))
        df$IsCell <- obs_df$star_IsCell

    #filter to at least some minimum number of transcripts.
    idx <- which(df$num_transcripts>=min_transcripts)
    df <- df[idx,]

    # filter sparse matrix to same barcodes in same order.
    dge <- sparse_matrix[,df$cell_barcode]

    result <- list(dge=dge, cell_features=df)
    return (result)
}

getOptimusGeneSymbols<-function (h5ad_file) {
    # Get the gene names.
    # Gene names can be repeated.
    # Read gene indices and corresponding gene names
    gene_indices <- h5read(h5ad_file, "/var/Gene")                      # Gene indices
    gene_names <- h5read(h5ad_file, "/var/__categories/Gene")           # Actual gene names

    # Map indices to names (if necessary)
    # If the indices are numeric, they map directly to the gene names
    gene_names_full <- gene_names[gene_indices + 1]
    return (as.character(gene_names_full))
}










##################
# JUNK
##################

# test_code_Optimus<-function () {
#     # Path to your h5ad file
#     h5ad_file <- "/broad/bican_um1_mccarroll/RNAseq/data/test/2023-08-25_v3_Bamboo_18d_10X_RNAseq_Optiprep8000_CaH_17k_rxn1/development_branch/2023-08-25_v3_Bamboo_18d_10X_RNAseq_Optiprep8000_CaH_17k_rxn1_dev_aug07_2024.h5ad"
#     h5ls(h5ad_file)
#     r=parseH5ad(h5ad_file)
#     r=prepareOptimusDataForNucleiSelection(r$dge, r$cell_features)
#     r$dge[1:5,1:5]
#
#     obs_df=r$cell_features
#     #filter to at least 20 transcripts
#     idx=which(obs_df$num_transcripts>=20)
#     obs_df=obs_df[idx,]
#     smoothScatter (log10(obs_df$num_transcripts), obs_df$pct_intronic, xlab="num_transcripts", ylab="% intronic", main="num_transcripts vs % intronic")
# }
#
# test_code_Optimus2<-function () {
#     library(rhdf5)
#     library(Matrix)
#
#     h5ad_file="/broad/mccarroll/nemesh/intronic_svm_cell_calling/svm_density/external_data/8k_pbmc.h5ad"
#     #h5ls(h5ad_file)
#     r=parseH5ad(h5ad_file)
#     r=prepareOptimusDataForNucleiSelection(r$dge, r$cell_features)
#
# }
#
# test_code_Optimus3<-function () {
#     library(rhdf5)
#     library(Matrix)
#     library (ggplot2)
#
#     h5ad_file="/broad/mccarroll/nemesh/intronic_svm_cell_calling/svm_density/external_data/macaque_test.h5ad"
#     #h5ls(h5ad_file)
#     #head (h5read(h5ad_file, cell_id_path))
#     r=parseH5ad(h5ad_file)
#     r=prepareOptimusDataForNucleiSelection(r$dge, r$cell_features)
#
#
#     p <- ggplot(r$cell_features, aes(x = log10(num_transcripts), y = pct_intronic, color = IsCell)) +
#         geom_point(size = 0.5, alpha=0.5) +
#         labs(
#             x = "log10(UMI)",
#             y = "% Intronic",
#             color = "STAR Selected Cell"
#         ) +
#         scale_color_manual(values = c("TRUE" = "green", "FALSE" = "lightblue")) +
#         theme_minimal() +
#         ggtitle(paste("STAR number of nuclei selected [", length(which(r$cell_features$IsCell==T)), "]")) +
#         guides(color = guide_legend(override.aes = list(size = 4, alpha=1)))
#
#
#     smoothScatter(log10(r$cell_features$num_transcripts), r$cell_features$pct_intronic, xlab="log10 UMIs", ylab="fraction intronic")
#
#
# }
#
#
#
# testParseCBRB<-function () {
#     h5ad_file <- "/broad/bican_um1_mccarroll/RNAseq/data/test/2024-06-12_v20_10X-GEX-3P_Pu_rxn1_optimus/submission_65efa5ef-3f35-4234-adea-7dc7594f1244/GRCh38_ensembl_v43.isa.exonic+intronic/cbrb/auto/v20_10X-GEX-3P_Pu_rxn1.h5"
#     #expression_matrix_path="/matrix"; gene_id_path="/matrix/features/name"; cell_id_path="/matrix/barcodes"
#     #file_structure <- h5ls(h5ad_file)
#     r=parseH5ad(h5ad_file, expression_matrix_path="/matrix", gene_id_path="/matrix/features/name", cell_id_path="/matrix/barcodes")
#     r$dge[1:5,1:5]
# }
