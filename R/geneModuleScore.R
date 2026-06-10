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
#' @param min_pseudobulk_observations The minimum number of pseudobulked
#'   observations required in each comparison group before differential
#'   expression and module scoring are attempted.
#' @param verbose A boolean indicating whether to print log messages.
#' @return A list containing the cell features with the gene module scores and
#'   QC plots. If there are no differentially expressed genes, the
#'   empty_gene_module_score will be set to NA and the plots will be NULL.
#' @import Seurat Matrix methods
#' @noRd
#'
computeSvmGeneModuleScore <- function(
  cell_features_labeled, dgeMatrix, numGenes = 100,
  useCellBenderFeatures = TRUE, negative_class = "empty",
  min_nucleus_exemplars = 10,
  min_negative_exemplars = 10,
  min_pseudobulk_observations = 10,
  verbose = FALSE
) {
  module_score_name <- paste0(negative_class, "_gene_module_score")

  log_info(
    "Computing gene module score [", module_score_name,
    "] for nucleus vs [", negative_class, "]"
  )

  empty_result <- makeEmptyGeneModuleResult(
    cell_features_labeled = cell_features_labeled,
    module_score_name = module_score_name
  )

  if (!"training_label_class" %in% colnames(cell_features_labeled)) {
    log_warn("training_label_class not found. Skipping gene module score.")
    return(empty_result)
  }

  num_nuclei <- sum(
    cell_features_labeled$training_label_class == "nucleus",
    na.rm = TRUE
  )
  num_negative <- sum(
    cell_features_labeled$training_label_class == negative_class,
    na.rm = TRUE
  )

  if (num_nuclei < min_nucleus_exemplars ||
    num_negative < min_negative_exemplars) {
    log_warn(
      "Skipping gene module score [", module_score_name,
      "] because class counts are nucleus [", num_nuclei,
      "] and ", negative_class, " [", num_negative,
      "]. Minimums are nucleus [", min_nucleus_exemplars,
      "] and ", negative_class, " [", min_negative_exemplars, "]."
    )
    return(empty_result)
  }

  dgeMatrix <- methods::as(dgeMatrix, "CsparseMatrix")
  gc()

  dgeMatrix <- dgeMatrix[, rownames(cell_features_labeled)]

  seurat_object <- CreateSeuratObject(counts = dgeMatrix)
  seurat_object <- add_cell_metadata(
    seurat_object = seurat_object,
    cell_features = cell_features_labeled,
    negative_class = negative_class
  )

  training_idx <- which(
    cell_features_labeled$training_label_class %in%
      c("nucleus", negative_class)
  )
  seurat_object_training <- seurat_object[, training_idx]

  seurat_object_pseudobulked <- pseudobulkNonTarget(
    seurat_object = seurat_object_training,
    negative_identity = negative_class,
    showPlot = FALSE,
    verbose = verbose
  )

  pseudobulk_group_counts <- countPseudobulkTrainingIdentities(
    seurat_object_pseudobulked
  )
  num_pseudobulk_nuclei <- getNamedCount(
    pseudobulk_group_counts,
    "nuclei"
  )
  num_pseudobulk_negative <- getNamedCount(
    pseudobulk_group_counts,
    negative_class
  )

  if (num_pseudobulk_nuclei < min_pseudobulk_observations ||
    num_pseudobulk_negative < min_pseudobulk_observations) {
    log_warn(
      "Skipping gene module score [", module_score_name,
      "] because pseudobulked class counts are nuclei [",
      num_pseudobulk_nuclei, "] and ", negative_class, " [",
      num_pseudobulk_negative, "]. Minimum required for each class is [",
      min_pseudobulk_observations, "]."
    )
    return(empty_result)
  }

  seurat_object_pseudobulked <- NormalizeData(seurat_object_pseudobulked,
    normalization.method = "LogNormalize", scale.factor = 10000,
    verbose = verbose
  )

  deWilcoxPB <- de_wilcox(
    seurat_object = seurat_object_pseudobulked,
    fdrThreshold = 0.05,
    min.pct = 0.25,
    logfc.threshold = 1,
    negative_identity = negative_class
  )
  geneListDown <- selectDifferentialGenes(deWilcoxPB, numGenes)

  if (checkNoDifferentialGenes(geneListDown)) {
    log_warn(
      "No differentially expressed genes found for module score [",
      module_score_name, "]."
    )
    return(empty_result)
  }

  moduleGeneTable <- makeModuleGeneTable(deWilcoxPB, geneListDown)

  seurat_object_pseudobulked <- AddModuleScore(seurat_object_pseudobulked,
    features = list(geneListDown), name = module_score_name,
    verbose = verbose
  )

  seurat_object <- NormalizeData(seurat_object,
    normalization.method = "LogNormalize", scale.factor = 10000,
    verbose = verbose
  )

  seurat_object <- AddModuleScore(seurat_object,
    features = list(geneListDown), name = module_score_name,
    verbose = verbose
  )

  scored_features <- extractCellFeatures(seurat_object)
  gene_module_plots <- generateGeneModulePlots(
    seurat_object_pseudobulked = seurat_object_pseudobulked,
    cell_features = scored_features,
    useCellBenderFeatures = useCellBenderFeatures,
    module_score_name = module_score_name,
    negative_class = negative_class
  )

  list(
    score = scored_features[[module_score_name]],
    downGenes = geneListDown,
    moduleGeneTable = moduleGeneTable,
    plots = gene_module_plots
  )
}

selectDifferentialGenes <- function(deWilcoxPB, numGenes) {
  downGenes <- deWilcoxPB[deWilcoxPB$avg_log2FC < 0, ]
  downGenes <- downGenes[order(downGenes$avg_log2FC, decreasing = FALSE), ]
  return(head(rownames(downGenes), numGenes))
}

makeModuleGeneTable <- function(deWilcoxPB, geneList) {
  if (length(geneList) == 0) {
    return(data.frame())
  }

  moduleGeneTable <- deWilcoxPB[geneList, , drop = FALSE]
  moduleGeneTable$gene <- rownames(moduleGeneTable)
  moduleGeneTable$rank <- seq_len(nrow(moduleGeneTable))
  rownames(moduleGeneTable) <- NULL

  moduleGeneTable[
    ,
    c("rank", "gene", setdiff(colnames(moduleGeneTable), c("rank", "gene"))),
    drop = FALSE
  ]
}

checkNoDifferentialGenes <- function(geneListDown) {
  return(length(geneListDown) == 0)
}

countPseudobulkTrainingIdentities <- function(seurat_object) {
  metadata <- seurat_object[[]]

  if (!"training_identity" %in% colnames(metadata)) {
    return(integer(0))
  }

  return(table(metadata$training_identity))
}

getNamedCount <- function(counts, name) {
  if (!name %in% names(counts)) {
    return(0L)
  }

  return(as.integer(counts[[name]]))
}

makeEmptyGeneModuleResult <- function(
  cell_features_labeled,
  module_score_name = "empty_gene_module_score"
) {
  list(
    score = rep(NA_real_, nrow(cell_features_labeled)),
    downGenes = NULL,
    moduleGeneTable = data.frame(),
    plots = makeEmptyGeneModulePlots(
      module_score_name = module_score_name
    )
  )
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
  cell_features, useCellBenderFeatures,
  module_score_name = "empty_gene_module_score",
  negative_class = "empty"
) {
  module_score_col <- paste0(module_score_name, "1")
  score_title <- paste0(negative_class, " score")

  training_plot <- VlnPlot(seurat_object_pseudobulked,
    features = c(module_score_col),
    group.by = "training_identity"
  )

  module_plot <- scatterPlotModuleScore(cell_features,
    moduleName = module_score_name,
    strTitle = score_title
  )

  result <- list()
  result[[paste0(module_score_name, "_training_data")]] <- training_plot
  result[[module_score_name]] <- module_plot

  if (useCellBenderFeatures && module_score_name == "empty_gene_module_score") {
    result[["frac_contamination"]] <- scatterPlotModuleScore(cell_features,
      moduleName = "frac_contamination",
      strTitle = "CellBender fraction contamination"
    )
    result[["empty_gene_module_score_vs_contam"]] <-
      plotModuleVsFracContam(cell_features,
        moduleName = module_score_name
      )
  }

  return(result)
}

makeEmptyGeneModulePlots <- function(
  module_score_name = "empty_gene_module_score",
  negative_class = "empty"
) {
  p_empty <- ggplot() +
    theme_void()

  result <- list()
  result[[paste0(module_score_name, "_training_data")]] <- p_empty
  result[[module_score_name]] <- p_empty

  if (module_score_name == "empty_gene_module_score") {
    result[["frac_contamination"]] <- p_empty
    result[["empty_gene_module_score_vs_contam"]] <- p_empty
  }

  result
}

de_wilcox <- function(
  seurat_object, fdrThreshold = 1, min.pct = 0.25,
  logfc.threshold = 0.25, negative_identity = "empty"
) {
  de_results <- FindMarkers(seurat_object,
    ident.1 = "nuclei",
    ident.2 = negative_identity, test.use = "wilcox", min.pct = min.pct,
    logfc.threshold = logfc.threshold, group.by = "training_identity",
    only.pos = FALSE
  )

  de_results <-
    de_results[order(de_results$avg_log2FC, decreasing = TRUE), ]

  de_results <- de_results[de_results$p_val_adj <= fdrThreshold, ]

  return(de_results)
}

add_cell_metadata <- function(
  seurat_object,
  cell_features,
  negative_class = "empty"
) {
  cell_features$training_identity <- NA_character_
  cell_features$training_identity[
    cell_features$training_label_class == "nucleus"
  ] <- "nuclei"
  cell_features$training_identity[
    cell_features$training_label_class == negative_class
  ] <- negative_class

  seurat_object <- AddMetaData(seurat_object, metadata = cell_features)

  Idents(seurat_object) <- cell_features$training_identity

  return(seurat_object)
}

#################
# PSEUDOBULKING
#################

#' Collapse empty droplets into larged pseudobulked droplets to match nuclei UMI
#' counts
#' @noRd
pseudobulkNonTarget <- function(
  seurat_object,
  negative_identity = "empty",
  showPlot = FALSE,
  verbose = FALSE
) {
  training_identity <- NULL

  nuclei_cells <- subset(seurat_object,
    subset = training_identity == "nuclei"
  )
  negative_cells <- subset(seurat_object,
    subset = training_identity == negative_identity
  )

  nuclei_expr <- GetAssayData(
    nuclei_cells,
    assay = "RNA",
    layer = "counts"
  )
  negative_expr <- GetAssayData(
    negative_cells,
    assay = "RNA",
    layer = "counts"
  )

  pseudobulked_negative_expr <-
    pseudobulk_to_match_umi_fast_optimized(nuclei_expr, negative_expr)

  pseudobulked_negative_expr <-
    removeZeroCountCells(pseudobulked_negative_expr)

  colnames(pseudobulked_negative_expr) <-
    paste0(toupper(negative_identity), "_", seq_len(
      ncol(pseudobulked_negative_expr)
    ))

  combined_counts <- cbind(nuclei_expr, pseudobulked_negative_expr)
  combined_seurat <- CreateSeuratObject(counts = combined_counts)

  combined_metadata <- createPseudobulkMetadata(
    nuclei_cells = nuclei_cells,
    pseudobulked_expr = pseudobulked_negative_expr,
    negative_identity = negative_identity
  )
  combined_seurat <- AddMetaData(combined_seurat, metadata = combined_metadata)

  if (showPlot) {
    plotPseudobulkDensity(combined_seurat)
  }

  log_info(
    "Finished constructing pseudobulked [", negative_identity,
    "] droplets"
  )
  return(combined_seurat)
}

removeZeroCountCells <- function(expression_matrix) {
  non_zero_cols <- colSums(expression_matrix) > 0
  return(expression_matrix[, non_zero_cols, drop = FALSE])
}

createPseudobulkMetadata <- function(
  nuclei_cells,
  pseudobulked_expr,
  negative_identity = "empty"
) {
  pseudobulked_cell_names <- colnames(pseudobulked_expr)

  pseudobulk_nCount_RNA <- colSums(pseudobulked_expr)
  pseudobulk_nFeature_RNA <- colSums(pseudobulked_expr > 0)
  pseudobulk_orig_ident <-
    rep("SeuratProject", length(pseudobulked_cell_names))

  empty_metadata <- nuclei_cells@meta.data[0, ]
  empty_metadata <- empty_metadata[rep(1, length(pseudobulked_cell_names)), ]

  rownames(empty_metadata) <- pseudobulked_cell_names
  empty_metadata$orig.ident <- pseudobulk_orig_ident
  empty_metadata$nCount_RNA <- pseudobulk_nCount_RNA
  empty_metadata$nFeature_RNA <- pseudobulk_nFeature_RNA
  empty_metadata$training_identity <- negative_identity

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
  non_target_expr, verbose = FALSE
) {
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
  cumsum_non_target, verbose
) {
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
  indices, non_target_indices, verbose
) {
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
