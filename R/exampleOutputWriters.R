#' Write Example SvmNucleusCallerInputs Data to 10x Format
#'
#' This function writes the `svmNucleusCallerInputs` dataset to a specified
#' temporary directory in 10x format. The output includes a compressed
#' Matrix Market (.mtx.gz) file, a features file (.tsv.gz) with gene metadata,
#' and a barcodes file (.tsv.gz) containing cell barcode information.
#'
#' @param output_dir The directory where the 10x-formatted files will be
#'   written.
#' @return A list of file paths for the generated files.
#' @examples
#' temp_dir <- tempdir()
#' writeExampleSvmNucleusCallerInputs(temp_dir)
#' list.files(temp_dir)  # Verify files were created
#' @export
writeExampleSvmNucleusCallerInputs <- function(output_dir) {
    # Load example dataset
    data(svmNucleusCallerInputs)

    # Ensure the directory exists
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    # Define file paths
    matrix_file <- file.path(output_dir, "matrix.mtx")
    features_file <- file.path(output_dir, "features.tsv")
    barcodes_file <- file.path(output_dir, "barcodes.tsv")

    # Convert DGE matrix to 10x format
    dge_matrix <- svmNucleusCallerInputs$dgeMatrix

    # Write Matrix Market formatted matrix
    Matrix::writeMM(dge_matrix, matrix_file)
    R.utils::gzip(matrix_file, overwrite = TRUE)

    # Write features (gene names) with 3 columns: gene symbol, gene symbol,
    # "Gene Expression"
    features <- data.frame(GENE = rownames(dge_matrix),
        SYMBOL = rownames(dge_matrix), TYPE = "Gene Expression")
    write.table(features, features_file, sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE)
    R.utils::gzip(features_file, overwrite = TRUE)

    # Write barcodes (cell barcodes)
    barcodes <- data.frame(CELL_BARCODE = colnames(dge_matrix))
    write.table(barcodes, barcodes_file, sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE)
    R.utils::gzip(barcodes_file, overwrite = TRUE)

    # Return the list of created files
    return(list(
        matrix_file = paste0(matrix_file, ".gz"),
        features_file = paste0(features_file, ".gz"),
        barcodes_file = paste0(barcodes_file, ".gz")
    ))
}

#' Write Example Cell Features to a File
#'
#' This function writes the `cellFeatures` from the `svmNucleusCallerInputs`
#' dataset to a specified directory in a tab-delimited format.
#'
#' @param output_dir The directory where the cell features file will be written.
#' @return The file path of the generated cell features file.
#' @examples
#' temp_dir <- tempdir()
#' writeExampleCellFeatures(temp_dir)
#' list.files(temp_dir)  # Verify file was created
#' @export
writeExampleCellFeatures <- function(output_dir) {
    # Load example dataset
    data(svmNucleusCallerInputs)

    # Ensure the directory exists
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    # Define file path
    cell_features_file <- file.path(output_dir, "cell_features.tsv")

    # Write cell features to file
    write.table(svmNucleusCallerInputs$cellFeatures, cell_features_file,
                sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    return(cell_features_file)
}

#' Write Example DGE Data as a Dense Flat File
#'
#' This function writes the `dgeMatrix` from the `svmNucleusCallerInputs`
#' dataset to a specified directory in a tab-delimited format. The first column
#' will be labeled "GENE" and contain gene names.
#'
#' @param output_dir The directory where the dense DGE file will be written.
#' @return The file path of the generated dense DGE file.
#' @examples
#' temp_dir <- tempdir()
#' writeExampleDenseDGE(temp_dir)
#' list.files(temp_dir)  # Verify file was created
#' @export
writeExampleDenseDGE <- function(output_dir) {
    # Load example dataset
    data(svmNucleusCallerInputs)

    # Ensure the directory exists
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    # Define file paths
    dge_file <- file.path(output_dir, "dge_matrix.tsv")
    dge_gz_file <- paste0(dge_file, ".gz")  # Path for the compressed file

    # Convert DGE matrix to data frame with "GENE" as the first column
    dge_matrix <- as.data.frame(as.matrix(svmNucleusCallerInputs$dgeMatrix))
    dge_matrix <- cbind(GENE = rownames(dge_matrix), dge_matrix)
    rownames(dge_matrix) <- NULL

    # Write dense DGE matrix to file
    write.table(dge_matrix, dge_file, sep = "\t", row.names = FALSE,
                col.names = TRUE, quote = FALSE)

    # Compress the file to .gz format
    R.utils::gzip(dge_file, overwrite = TRUE)

    # Return the final compressed file name
    return(dge_gz_file)
}
