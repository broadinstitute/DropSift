###################
# GENE MODULE AND FEATURE BUILDING BASED ON EXPRESSION
###################

#' Find Nuclei/Empty gene modules
#'
#' Given the training nuclei and empty cell barcodes, find the gene modules that
#' are more differentially expressed between the groups, and score all cells
#' with those features.
#'
#' @param cell_features_labeled The cell_features_labelel dataframe with an
#'   additional indicator column 'training_label_is_cell' that is TRUE for
#'   nuclei and FALSE for empty.  This value is set to NA for all other cell
#'   barcodes.
#' @param dgeMatrix The matrix of gene expression.
#' @param numGenes The number of genes to select for each module.
#' @param useCellBenderFeatures If true, additional plots are generated to
#'   compare the empty module score to cellbender fraction of UMIs removed.
#' @param verbose A boolean indicating whether to print log messages.
#' @return A list containing the cell features with the gene module scores and
#'   QC plots. If there are no differentially expressed genes, the
#'   empty_gene_module_score will be set to NA and the plots will be NULL.
#' @import Seurat Matrix methods
#' @noRd
#'
addGeneModules <- function(
    cell_features_labeled, dgeMatrix, numGenes = 100,
    useCellBenderFeatures = TRUE, verbose = FALSE) {
  log_info("Adding gene module score(s)")
  dgeMatrix <- methods::as(dgeMatrix, "CsparseMatrix")
  gc() # Garbage collection for memory efficiency
  # Ensure dgeMatrix is ordered correctly
  dgeMatrix <- dgeMatrix[, rownames(cell_features_labeled)]

  # Create a Seurat object with expression and metadata
  seurat_object <- CreateSeuratObject(counts = dgeMatrix)
  seurat_object <- add_cell_metadata(seurat_object, cell_features_labeled)

  # Subset training data only
  seurat_object_trainng <-
    seurat_object[, !is.na(seurat_object$training_label_is_cell)]

  # Pseudobulk empty droplets for differential expression
  seurat_object_pseudobulked <- pseudobulkEmpties(seurat_object_trainng,
    showPlot = FALSE, verbose = verbose
  )
  seurat_object_pseudobulked <- NormalizeData(seurat_object_pseudobulked,
    normalization.method = "LogNormalize", scale.factor = 10000,
    verbose = verbose
  )

  # Get differentially expressed genes
  deWilcoxPB <- de_wilcox(seurat_object_pseudobulked,
    fdrThreshold = 0.05,
    min.pct = 0.25, logfc.threshold = 1
  )
  geneListDown <- selectDifferentialGenes(deWilcoxPB, numGenes)

  if (checkNoDifferentialGenes(geneListDown)) {
    return(handleNoDifferentialGenes(cell_features_labeled))
  }

  # Score the module in pseudobulk and full dataset
  seurat_object_pseudobulked <- AddModuleScore(seurat_object_pseudobulked,
    features = list(geneListDown), name = "empty_gene_module_score",
    verbose = verbose
  )
  seurat_object <- NormalizeData(seurat_object,
    normalization.method = "LogNormalize", scale.factor = 10000,
    verbose = verbose
  )

  seurat_object <- AddModuleScore(seurat_object,
    features = list(geneListDown), name = "empty_gene_module_score",
    verbose = verbose
  )

  # Extract cell features and clean column names
  cell_features <- extractCellFeatures(seurat_object)
  gene_module_plots <- generateGeneModulePlots(
    seurat_object_pseudobulked,
    cell_features, useCellBenderFeatures
  )
  return(list(
    cell_features = cell_features, downGenes = geneListDown,
    gene_module_plots = gene_module_plots
  ))
}

selectDifferentialGenes <- function(deWilcoxPB, numGenes) {
  downGenes <- deWilcoxPB[deWilcoxPB$avg_log2FC < 0, ]
  downGenes <- downGenes[order(downGenes$avg_log2FC, decreasing = FALSE), ]
  return(head(rownames(downGenes), numGenes))
}

checkNoDifferentialGenes <- function(geneListDown) {
  return(length(geneListDown) == 0)
}

handleNoDifferentialGenes <- function(cell_features_labeled) {
  log_warn(
    "No differentially expressed genes found",
    "for empty cell barcodes."
  )

  return(list(
    cell_features = cell_features_labeled, downGenes = NULL,
    gene_module_plots = list(
      training_data = NULL,
      frac_contamination = NULL, empty_gene_module_score = NULL,
      empty_gene_module_score_vs_contam = NULL
    )
  ))
}

extractCellFeatures <- function(seurat_object) {
  cell_features <- seurat_object[[]]

  colnames(cell_features) <- gsub(
    x = colnames(cell_features),
    pattern = "1", replacement = ""
  )

  dropCols <- c(
    "orig.ident", "nCount_RNA", "nFeature_RNA",
    "training_identity"
  )

  return(cell_features[, !(colnames(cell_features) %in% dropCols)])
}

generateGeneModulePlots <- function(
    seurat_object_pseudobulked,
    cell_features, useCellBenderFeatures) {
  p1 <- VlnPlot(seurat_object_pseudobulked,
    features = c("empty_gene_module_score1"),
    group.by = "training_identity"
  )
  p2 <- ggplot() +
    theme_void() # Empty plot
  p6 <- ggplot() +
    theme_void() # Empty plot

  p4 <- scatterPlotModuleScore(cell_features,
    moduleName = "empty_gene_module_score",
    strTitle = "Empty Score"
  )

  if (useCellBenderFeatures) {
    p2 <- scatterPlotModuleScore(cell_features,
      moduleName = "frac_contamination",
      strTitle = "CellBender fraction contamination"
    )
    p6 <- plotModuleVsFracContam(cell_features,
      moduleName = "empty_gene_module_score"
    )
  }

  return(list(
    training_data = p1, frac_contamination = p2,
    empty_gene_module_score = p4, empty_gene_module_score_vs_contam = p6
  ))
}

de_wilcox <- function(
    seurat_object, fdrThreshold = 1, min.pct = 0.25,
    logfc.threshold = 0.25) {
  de_results <- FindMarkers(seurat_object,
    ident.1 = "nuclei",
    ident.2 = "empty", test.use = "wilcox", min.pct = min.pct,
    logfc.threshold = logfc.threshold, group.by = "training_identity",
    only.pos = FALSE
  )

  de_results <-
    de_results[order(de_results$avg_log2FC, decreasing = TRUE), ]

  de_results <- de_results[de_results$p_val_adj <= fdrThreshold, ]

  return(de_results)
}

add_cell_metadata <- function(seurat_object, cell_features) {
  cell_features$training_identity <-
    ifelse(cell_features$training_label_is_cell, "nuclei", "empty")

  seurat_object <- AddMetaData(seurat_object, metadata = cell_features)

  Idents(seurat_object) <- cell_features$training_label_is_cell

  return(seurat_object)
}

#################
# PSEUDOBULKING
#################

#' Collapse empty droplets into larged pseudobulked droplets to match nuclei UMI
#' counts
#' @noRd
pseudobulkEmpties <- function(
    seurat_object, showPlot = FALSE,
    verbose = FALSE) {
  training_identity <- NULL # for R CMD CHECK
  # Separate nuclei and empty droplets
  nuclei_cells <- subset(seurat_object,
    subset = training_identity == "nuclei"
  )
  empty_cells <- subset(seurat_object,
    subset = training_identity == "empty"
  )

  # Extract expression matrices
  nuclei_expr <- GetAssayData(nuclei_cells, assay = "RNA", layer = "counts")
  empty_expr <- GetAssayData(empty_cells, assay = "RNA", layer = "counts")

  # Pseudobulk empty droplets to match UMI counts of nuclei
  pseudobulked_empty_expr <-
    pseudobulk_to_match_umi_fast_optimized(nuclei_expr, empty_expr)

  # Remove columns with zero total counts
  pseudobulked_empty_expr <- removeZeroCountCells(pseudobulked_empty_expr)

  # Assign unique names to pseudobulked empty droplets
  colnames(pseudobulked_empty_expr) <-
    paste0("EMPTY_", seq_len(ncol(pseudobulked_empty_expr)))

  # Combine nuclei and pseudobulked empty data
  combined_counts <- cbind(nuclei_expr, pseudobulked_empty_expr)
  combined_seurat <- CreateSeuratObject(counts = combined_counts)

  # Generate and merge metadata
  combined_metadata <-
    createPseudobulkMetadata(nuclei_cells, pseudobulked_empty_expr)
  combined_seurat <-
    AddMetaData(combined_seurat, metadata = combined_metadata)

  # Generate QC plots if required
  if (showPlot) plotPseudobulkDensity(combined_seurat)

  log_info("Finished constructing pseudobulked empty droplets")
  return(combined_seurat)
}

removeZeroCountCells <- function(expression_matrix) {
  non_zero_cols <- colSums(expression_matrix) > 0
  return(expression_matrix[, non_zero_cols, drop = FALSE])
}

createPseudobulkMetadata <- function(nuclei_cells, pseudobulked_empty_expr) {
  pseudobulked_cell_names <- colnames(pseudobulked_empty_expr)

  # Compute metadata for pseudobulked cells
  pseudobulk_nCount_RNA <- colSums(pseudobulked_empty_expr)
  pseudobulk_nFeature_RNA <- colSums(pseudobulked_empty_expr > 0)
  pseudobulk_orig_ident <-
    rep("SeuratProject", length(pseudobulked_cell_names))

  # Create an empty metadata frame matching nuclei_cells@meta.data
  empty_metadata <- nuclei_cells@meta.data[0, ]
  empty_metadata <- empty_metadata[rep(1, length(pseudobulked_cell_names)), ]

  # Assign row names and fill in metadata fields
  rownames(empty_metadata) <- pseudobulked_cell_names
  empty_metadata$orig.ident <- pseudobulk_orig_ident
  empty_metadata$nCount_RNA <- pseudobulk_nCount_RNA
  empty_metadata$nFeature_RNA <- pseudobulk_nFeature_RNA
  empty_metadata$training_identity <- "empty"

  return(rbind(nuclei_cells@meta.data, empty_metadata))
}

plotPseudobulkDensity <- function(seurat_object) {
  z <- seurat_object@meta.data

  # Ensure variables are explicitly defined to avoid R CMD CHECK issues
  nCount_RNA <- training_identity <- NULL

  p <- ggplot(z, aes(x = nCount_RNA, fill = training_identity)) +
    geom_density(alpha = 0.5) +
    labs(
      title = "Density Plot of nCount_RNA by Training Identity",
      x = "nCount_RNA", y = "Density"
    ) +
    theme_minimal()

  gridExtra::grid.arrange(p)
}


pseudobulk_to_match_umi_fast_optimized <- function(
    target_expr,
    non_target_expr, verbose = FALSE) {
  log_info(
    "Constructing pseudobulked empty ",
    "droplets for gene module detection"
  )

  # Compute UMI counts for nuclei and empty droplets
  target_umis <- Matrix::colSums(target_expr)
  shuffled_target_umis <- sample(target_umis) # Shuffle target UMIs

  non_target_umi_counts <- Matrix::colSums(non_target_expr)
  non_target_indices <- sample(seq_along(non_target_umi_counts))
  shuffled_non_target_umis <- non_target_umi_counts[non_target_indices]

  # Initialize sparse pseudobulk matrix
  pseudobulk_matrix <- initializePseudobulkMatrix(
    non_target_expr,
    length(shuffled_target_umis)
  )

  # Compute cumulative sums
  cumsum_target_umis <- cumsum(shuffled_target_umis)
  cumsum_non_target <- cumsum(shuffled_non_target_umis)

  # Determine start and end indices for pseudobulking
  indices <- computePseudobulkIndices(
    cumsum_target_umis,
    cumsum_non_target, verbose
  )

  # Aggregate non-target expression to match target UMIs
  pseudobulk_matrix <- aggregatePseudobulkedUMIs(
    pseudobulk_matrix,
    non_target_expr, indices, non_target_indices, verbose
  )

  return(pseudobulk_matrix)
}

initializePseudobulkMatrix <- function(non_target_expr, num_target_cells) {
  pseudobulk_matrix <- Matrix::Matrix(0,
    nrow = nrow(non_target_expr),
    ncol = num_target_cells, sparse = TRUE
  )

  rownames(pseudobulk_matrix) <- rownames(non_target_expr)
  return(pseudobulk_matrix)
}

computePseudobulkIndices <- function(
    cumsum_target_umis,
    cumsum_non_target, verbose) {
  start_indices <- integer(length(cumsum_target_umis))
  end_indices <- integer(length(cumsum_target_umis))
  prev_end_index <- 0

  for (i in seq_along(cumsum_target_umis)) {
    end_index <- findInterval(cumsum_target_umis[i], cumsum_non_target) + 1

    if (end_index > length(cumsum_non_target)) {
      if (verbose) {
        message(
          "Not enough non-target UMIs",
          "to match target UMIs."
        )
      }
      break
    }

    start_indices[i] <- prev_end_index + 1
    end_indices[i] <- end_index
    prev_end_index <- end_index
  }

  valid <- which(start_indices > 0 & end_indices >= start_indices)
  return(list(start = start_indices[valid], end = end_indices[valid]))
}

aggregatePseudobulkedUMIs <- function(
    pseudobulk_matrix, non_target_expr,
    indices, non_target_indices, verbose) {
  for (i in seq_along(indices$start)) {
    selected_cells_indices <-
      non_target_indices[indices$start[i]:indices$end[i]]

    pseudobulk_matrix[, i] <-
      Matrix::rowSums(
        non_target_expr[, selected_cells_indices, drop = FALSE]
      )

    if (verbose && i %% 100 == 0) {
      count <- ncol(pseudobulk_matrix)
      message("Pseudobulked cell barcode [", i, "] of [", count, "]")
    }
  }
  return(pseudobulk_matrix)
}
