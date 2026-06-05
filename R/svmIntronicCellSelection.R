# Some global parameters
MIN_UMIs_PER_STAMP <- 20 # what's worth even including in analysis
log10_UMI_AXIS_RANGE_NEW <- c(log10(MIN_UMIs_PER_STAMP), 6) # for plotting

#' Call nuclei using an SVM trained on cell summary features.
#'
#' A convenience method for calling cells using an SVM trained on cell summary
#' features.
#'
#' Please set a random seed for reproducibility, e.g., `set.seed(1)` or by
#' passing a value to random.seed parameter. This ensures that the SVM
#' initialization and training are consistent across runs.
#'
#' @param datasetName The name of the dataset.  Used for plotting.  Useful to
#'   use a unique experimental identifier (UEI) if available.
#' @param cellFeaturesFile The cell features file.  This file must contain a
#'   column cell_barcode that is used to match the cell barcodes in the DGE
#'   matrix and a num_reads column that defines the number of uniquely mapped
#'   reads for the cell.  Additional expected columns used to train the SVM are
#'   defined by the features vector.
#' @param dgeMatrixFile A file containing a dense DGE matrix in DropSeq format,
#'   or the directory where files in the MTX format exist. If a directory is
#'   supplied, 3 files are expect in the directory: matrix.mtx.gz,
#'   features.tsv.gz, barcodes.tsv.gz that encode the counts in sparse format,
#'   the genes, and the cell barcodes respectively.
#' @param optimusH5File The h5ad file produced by the Optimus pipeline.  This
#'   input contains both the DGE matrix and the cell features required for cell
#'   calling.  Cell level metrics are appropriately transformed to match SVM
#'   expectations for training.
#' @param cellProbabilityThreshold The probability threshold for selecting
#'   cells.  This value scales from 0-1, with 1 being extremely
#'   confident/stringent.  Set to null to use the SVM defaults (>0.5 = nuclei).
#'   This modifies plotting outputs and the output features file classification,
#'   but does not filter any cells from the output.
#' @param max_umis_empty All cell barcodes with fewer than this number of UMIs
#'   are not considered for any aspects of nuclei selection.  It is assumed
#'   these cell barcodes have fewer UMIs than the ambient RNA.
#' @param features A list of features to use for cell selection.  The features
#'   are used to train the SVM.  These features are columns of the cell features
#'   data frame, with the exception of empty_gene_module_score, which is
#'   generated from the gene expression data.
#' @param useCBRBFeatures When true, the cell bender feature frac_contamination
#'   is used for cell selection. When false, these features are not used.  This
#'   modifies the features argument.
#' @param useCBRBInitialization When true, the cellbender feature
#' frac_contamination is used to select exemplar nuclei and empty droplets.
#' This option can be false and useCBRBFeatures to force DropSift to use the
#' non-CBRB initialization but still include CBRB in the model features.
#' @param forceTwoClusterSolution When true, the initialization of the SVM will
#'   attempt to find a solution with two clusters. In cases where an experiment
#'   is overloaded with nuclei, this may correct the initial set of nuclei and
#'   empty droplets selected. This argument is specific to the density-style
#'   selection of exemplars that is only applicable when useCBRBFeatures is
#'   false.
#' @param use2DTrainingRefinement EXPERIMENTAL (!!!) When true, the density-only
#'   initialization refines rectangular empty-droplet and nucleus exemplar
#'   selections with connected components from a two-dimensional HDR. The
#'   default is false.
#' @param outPDF The PDF file to write the plots to.
#' @param outFeaturesFile The cell features dataframe, further annotated by the
#'   SVM to include the cell probability and label, along with which cell
#'   barcodes were used for training.
#' @param outCellBenderInitialParameters If useCBRBFeatures is false and this
#'   parameter is not null, the SVM is trained without cellbenber features, and
#'   the parameters expected_cells and total_droplets_included are
#' estimated and written to this file.
#' @param random.seed The random seed to use for reproducibility.  Default is NA, which means no seed is set.
#' @return This function does not return a value. Instead, the results are
#'   written to the specified output files (`outPDF`, `outFeaturesFile`, and
#'   `outCellBenderInitialParameters`).
#'
#' @import logger
#' @importFrom utils write.table
#' @export
#' @examples
#' # Create a temporary directory for test data
#' temp_dir <- tempdir()
#'
#' # Write example data to files
#' dge_files <- writeExampleSvmNucleusCallerInputs(temp_dir)
#' cell_features_file <- writeExampleCellFeatures(temp_dir)
#'
#' # Define output file paths
#' out_pdf_file <- file.path(temp_dir, "output.pdf")
#' out_features_file <- file.path(temp_dir, "classified_features.tsv")
#'
#' # Run the function
#' runIntronicSVM(
#'   datasetName = "example_dataset",
#'   cellFeaturesFile = cell_features_file,
#'   dgeMatrixFile = temp_dir, # directory with compressed 10x files
#'   useCBRBFeatures = FALSE,
#'   forceTwoClusterSolution = FALSE,
#'   outPDF = out_pdf_file,
#'   outFeaturesFile = out_features_file,
#'   random.seed = 1
#' )
#'
#' # Check that output files were generated
#' file.exists(out_pdf_file) # Should return TRUE
#' file.exists(out_features_file) # Should return TRUE
#'
#' # Inspect the PDF output report
#' if (file.exists(out_pdf_file)) {
#'   utils::browseURL(out_pdf_file)
#' }
runIntronicSVM <- function(
  datasetName, cellFeaturesFile = NULL,
  dgeMatrixFile = NULL, optimusH5File = NULL, cellProbabilityThreshold = NULL,
  max_umis_empty = 50, features = NULL, useCBRBFeatures = TRUE,
  useCBRBInitialization = useCBRBFeatures,
  forceTwoClusterSolution = FALSE, use2DTrainingRefinement = FALSE,
  outPDF = NULL, outFeaturesFile = NULL,
  outCellBenderInitialParameters = NULL,
  random.seed = NA
) {
  if (!is.na(random.seed)) {
    set.seed(random.seed)
  }

  # parse and validate inputs.
  r <- parseInputs(
    cellFeaturesFile = cellFeaturesFile,
    dgeMatrixFile = dgeMatrixFile,
    optimusH5File = optimusH5File
  )

  cell_features <- r$cell_features
  dgeMatrix <- r$dgeMatrix
  svmNucleusCaller <- SvmNucleusCaller(
    cellFeatures = cell_features,
    dgeMatrix = dgeMatrix,
    cellProbabilityThreshold = cellProbabilityThreshold,
    maxUmisEmpty = max_umis_empty, featureColumns = features,
    forceTwoClusterSolution = forceTwoClusterSolution,
    use2DTrainingRefinement = use2DTrainingRefinement,
    useCBRBFeatures = useCBRBFeatures,
    useCBRBInitialization = useCBRBInitialization,
    datasetName = datasetName
  )

  cell_features_result <- svmNucleusCaller$cell_features

  # Open the PDF device
  if (!is.null(outPDF)) {
    grDevices::pdf(outPDF)
  }

  plotSvmNucleusCaller(svmNucleusCaller)

  if (!is.null(outPDF)) {
    grDevices::dev.off()
  }

  if (!is.null(outFeaturesFile)) {
    write.table(cell_features_result, outFeaturesFile,
      row.names = FALSE,
      col.names = TRUE, quote = FALSE, sep = "\t"
    )
  }

  # write initial cellbender remove-background parameters to file
  if (!is.null(outCellBenderInitialParameters)) {
    cbrbArgs <- getCBRBArgs(svmNucleusCaller)
    result <- data.frame(cbrbArgs)
    write.table(result, outCellBenderInitialParameters,
      row.names = FALSE,
      col.names = TRUE, quote = FALSE, sep = "\t"
    )
  }
}

#' Select cells using an SVM trained on cell summary features
#'
#' @param dataset_name The name of the dataset.  Used for plotting.  Useful to
#'   use a unique experimental identifier (UEI) if available.
#' @param cell_features A dataframe of cell features produced by
#'   buildCellFeaturesSimple()
#' @param dgeMatrix A dense or sparse DGE matrix.
#' @inheritParams runIntronicSVM
#' @return A list containing the dataset name, the cell features with the SVM
#'   results, various plots, and the DGE matrix if it was parsed.
#' @import hdrcde ggplot2 grid gridExtra cowplot ggrastr logger gridGraphics
#' @importFrom e1071 svm
#' @importFrom stats predict
#' @noRd
callByIntronicSVM <- function(
  dataset_name, cell_features, dgeMatrix,
  cellProbabilityThreshold = NULL, max_umis_empty = 50,
  features, useCBRBInitialization = TRUE, forceTwoClusterSolution = FALSE,
  use2DTrainingRefinement = FALSE
) {
  validateFeaturePresence(cell_features, features)
  maxContaminationThreshold <- 0.1 # CBRB-specific contamination threshold
  log_info("Beginning SVM Intronic Cell Selection")
  featureListStr <- paste(features, collapse = ", ")
  log_info("Features used for cell selection:", featureListStr)
  rownames(cell_features) <- cell_features$cell_barcode
  cell_features$cell_barcode <- NULL
  useCellBenderFeatures <- "frac_contamination" %in% features
  # don't use useCBRBInitialization without cellbender features
  if (!useCellBenderFeatures) {
    useCBRBInitialization <- FALSE
  }
  allBounds <- findTrainingDataBounds(cell_features, max_umis_empty,
    useCellBenderFeatures = useCBRBInitialization,
    forceTwoClusterSolution = forceTwoClusterSolution,
    use2DTrainingRefinement = use2DTrainingRefinement
  )
  bounds_empty <- allBounds$bounds_empty
  bounds_non_empty <- allBounds$bounds_non_empty
  cell_features_labeled <- labelTrainingData(
    cell_features = cell_features,
    bounds_empty = bounds_empty,
    bounds_non_empty = bounds_non_empty,
    maxContaminationThreshold = maxContaminationThreshold,
    useCellBenderFeatures = useCellBenderFeatures,
    training_empty_barcodes = allBounds$training_empty_barcodes,
    training_nucleus_barcodes = allBounds$training_nucleus_barcodes,
    training_debris_barcodes = allBounds$training_debris_barcodes,
    bounds_debris = allBounds$bounds_debris
  )
  logTrainingDataSelection(cell_features_labeled)

  empty_module_result <- computeSvmGeneModuleScore(
    cell_features_labeled = cell_features_labeled,
    dgeMatrix = dgeMatrix,
    numGenes = 100,
    useCellBenderFeatures = useCellBenderFeatures,
    negative_class = "empty",
    min_nucleus_exemplars = 100,
    min_negative_exemplars = 100,
    verbose = FALSE
  )

  cell_features_labeled$empty_gene_module_score <-
    empty_module_result$score

  debris_module_result <- computeSvmGeneModuleScore(
    cell_features_labeled = cell_features_labeled,
    dgeMatrix = dgeMatrix,
    numGenes = 100,
    useCellBenderFeatures = useCellBenderFeatures,
    negative_class = "debris",
    min_nucleus_exemplars = 100,
    min_negative_exemplars = 100,
    verbose = FALSE
  )

  cell_features_labeled$debris_gene_module_score <-
    debris_module_result$score

  geneModulePlots <- c(
    empty_module_result$plots,
    debris_module_result$plots
  )

  svm_empty_result <- runBinarySVM(
    cell_features_labeled = cell_features_labeled,
    features = features,
    positive_class = "nucleus",
    negative_class = "empty",
    probability_col = "p_nucleus_vs_empty",
    min_positive_exemplars = 100,
    min_negative_exemplars = 100
  )

  if (is.null(svm_empty_result)) {
    stop("Unable to train the nucleus-vs-empty SVM.")
  }

  features_nucleus_vs_debris <- c(
    setdiff(features, "empty_gene_module_score"),
    "debris_gene_module_score"
  )

  svm_debris_result <- NULL
  if (!is.null(debris_module_result$downGenes)) {
    svm_debris_result <- runBinarySVM(
      cell_features_labeled = cell_features_labeled,
      features = features_nucleus_vs_debris,
      positive_class = "nucleus",
      negative_class = "debris",
      probability_col = "p_nucleus_vs_debris",
      min_positive_exemplars = 100,
      min_negative_exemplars = 100
    )
  }

  features_debris_vs_empty <- unique(c(
    features,
    "debris_gene_module_score"
  ))

  svm_debris_vs_empty_result <- NULL
  if (!is.null(debris_module_result$downGenes)) {
    svm_debris_vs_empty_result <- runBinarySVM(
      cell_features_labeled = cell_features_labeled,
      features = features_debris_vs_empty,
      positive_class = "debris",
      negative_class = "empty",
      probability_col = "p_debris_vs_empty",
      min_positive_exemplars = 100,
      min_negative_exemplars = 100
    )
  }

  cell_features_result <- cell_features_labeled
  cell_features_result$p_nucleus_vs_empty <-
    svm_empty_result$probabilities

  gene_module_exemplar_plot <- NULL

  if (is.null(svm_debris_result)) {
    cell_features_result$p_nucleus_vs_debris <- NA_real_
    cell_features_result$is_cell_prob <-
      cell_features_result$p_nucleus_vs_empty
    cell_features_result$p_debris_vs_nucleus <- NA_real_
  } else {
    gene_module_exemplar_plot <- plotGeneModuleScoresByExemplarClass(
      cell_features_labeled,
      strTitle = "Gene module scores by exemplar class"
    )

    cell_features_result$p_nucleus_vs_debris <-
      svm_debris_result$probabilities

    cell_features_result$is_cell_prob <- pmin(
      cell_features_result$p_nucleus_vs_empty,
      cell_features_result$p_nucleus_vs_debris
    )

    cell_features_result$p_debris_vs_nucleus <-
      1 - cell_features_result$p_nucleus_vs_debris
  }

  if (is.null(svm_debris_vs_empty_result)) {
    cell_features_result$p_debris_vs_empty <- NA_real_
    cell_features_result$is_debris_prob <- NA_real_
  } else {
    cell_features_result$p_debris_vs_empty <-
      svm_debris_vs_empty_result$probabilities

    cell_features_result$is_debris_prob <- pmin(
      cell_features_result$p_debris_vs_empty,
      cell_features_result$p_debris_vs_nucleus
    )
  }

  probabilityThreshold <- 0.5
  if (!is.null(cellProbabilityThreshold)) {
    probabilityThreshold <- cellProbabilityThreshold
  }

  cell_features_result$barcode_class <- assignBarcodeClass(
    is_cell_prob = cell_features_result$is_cell_prob,
    is_debris_prob = cell_features_result$is_debris_prob,
    probabilityThreshold = probabilityThreshold
  )

  # This is a post-SVM rule, not part of any binary classifier. Barcodes
  # below the upper UMI edge of the empty exemplar region are not assigned a
  # final barcode class, even if an SVM assigns them nonzero probability.
  idx_low_umi <- which(log10(cell_features_result$num_transcripts + 1) <
    bounds_empty$umi_upper_bound)

  cell_features_result$is_cell_prob[idx_low_umi] <- NA_real_
  cell_features_result$is_debris_prob[idx_low_umi] <- NA_real_
  cell_features_result$barcode_class[idx_low_umi] <- NA_character_

  selectionPlots <- createSelectionVisualization(
    cell_features_result,
    bounds_empty, bounds_non_empty, useCellBenderFeatures, dataset_name,
    bounds_debris = allBounds$bounds_debris,
    training_empty_barcodes = allBounds$training_empty_barcodes,
    training_nucleus_barcodes = allBounds$training_nucleus_barcodes,
    training_debris_barcodes = allBounds$training_debris_barcodes,
    umi_threshold = allBounds$best_umi_threshold,
    debris_intronic_floor = allBounds$debris_intronic_floor
  )
  feature_plot_nucleus_vs_empty <- plotScaledTrainingDataFeatures(
    svm_empty_result$trainingData,
    strTitle = "Training Features: nucleus vs empty"
  )
  feature_plot_nucleus_vs_debris <- NULL
  if (!is.null(svm_debris_result)) {
    feature_plot_nucleus_vs_debris <- plotScaledTrainingDataFeatures(
      svm_debris_result$trainingData,
      strTitle = "Training Features: nucleus vs debris"
    )
  }

  log_info("Nuclei selection finished")
  return(list(
    dataset_name = dataset_name,
    cell_features = cell_features_result,
    features = features,
    features_nucleus_vs_debris = features_nucleus_vs_debris,
    features_debris_vs_empty = features_debris_vs_empty,
    plots = selectionPlots,
    featurePlotNucleusVsEmpty = feature_plot_nucleus_vs_empty,
    featurePlotNucleusVsDebris = feature_plot_nucleus_vs_debris,
    geneModulePlots = geneModulePlots,
    geneModuleExemplarPlot = gene_module_exemplar_plot,
    empty_down_genes = empty_module_result$downGenes,
    debris_down_genes = debris_module_result$downGenes,
    empty_gene_module_genes = empty_module_result$moduleGeneTable,
    debris_gene_module_genes = debris_module_result$moduleGeneTable,
    bounds_empty = bounds_empty, bounds_non_empty = bounds_non_empty,
    bounds_debris = allBounds$bounds_debris,
    training_empty_barcodes = allBounds$training_empty_barcodes,
    training_nucleus_barcodes = allBounds$training_nucleus_barcodes,
    training_debris_barcodes = allBounds$training_debris_barcodes,
    svm_model_nucleus_vs_empty = svm_empty_result$svm_model,
    svm_model_nucleus_vs_debris = if (is.null(svm_debris_result)) {
      NULL
    } else {
      svm_debris_result$svm_model
    },
    svm_model_debris_vs_empty = if (is.null(svm_debris_vs_empty_result)) {
      NULL
    } else {
      svm_debris_vs_empty_result$svm_model
    },
    trainingData_nucleus_vs_empty = svm_empty_result$trainingData,
    trainingData_nucleus_vs_debris = if (is.null(svm_debris_result)) {
      NULL
    } else {
      svm_debris_result$trainingData
    },
    trainingData_debris_vs_empty = if (is.null(svm_debris_vs_empty_result)) {
      NULL
    } else {
      svm_debris_vs_empty_result$trainingData
    },
    use2DTrainingRefinement = use2DTrainingRefinement
  ))
}

assignBarcodeClass <- function(
  is_cell_prob,
  is_debris_prob,
  probabilityThreshold = 0.5
) {
  barcode_class <- rep("empty_or_other", length(is_cell_prob))

  idxNucleus <- !is.na(is_cell_prob) & is_cell_prob >= probabilityThreshold
  idxDebris <- !is.na(is_debris_prob) &
    is_debris_prob >= probabilityThreshold

  barcode_class[idxNucleus] <- "nucleus"
  barcode_class[idxDebris] <- "debris"

  idxBoth <- which(idxNucleus & idxDebris)
  if (length(idxBoth) > 0) {
    barcode_class[idxBoth] <- ifelse(
      is_cell_prob[idxBoth] >= is_debris_prob[idxBoth],
      "nucleus",
      "debris"
    )
  }

  barcode_class
}

validateFeaturePresence <- function(cell_features, features) {
  for (f in setdiff(features, "empty_gene_module_score")) {
    if (!(f %in% colnames(cell_features))) {
      stop(paste(
        "Requested Features [", paste(features, collapse = ", "),
        "]. Feature", f, "not found in cell_features data frame."
      ))
    }
  }
}

createSelectionVisualization <- function(
  cell_features_labeled, bounds_empty,
  bounds_non_empty, useCellBenderFeatures, dataset_name,
  bounds_debris = NULL, training_empty_barcodes = NULL,
  training_nucleus_barcodes = NULL, training_debris_barcodes = NULL,
  umi_threshold = NULL, debris_intronic_floor = NULL
) {
  p1 <- plotExpressionVsIntronic(cell_features_labeled,
    title = "All cell barcodes",
    useCellBenderFeatures = useCellBenderFeatures
  )

  p2 <- cowplot::ggdraw(function() {
    plotCellTypeIntervals(
      cell_features_labeled, bounds_empty,
      bounds_non_empty,
      bounds_debris = bounds_debris,
      training_empty_barcodes = training_empty_barcodes,
      training_nucleus_barcodes = training_nucleus_barcodes,
      training_debris_barcodes = training_debris_barcodes,
      umi_threshold = umi_threshold,
      debris_intronic_floor = debris_intronic_floor,
      show_1d_bounds = is.null(training_empty_barcodes) ||
        is.null(training_nucleus_barcodes),
      show_2d_bounds = !is.null(training_empty_barcodes) &&
        !is.null(training_nucleus_barcodes)
    )
  })

  p3 <- plotSelectedCells(cell_features_labeled)
  p4 <- plotCellProbabilities(cell_features_labeled,
    strTitle = "Cell Probability"
  )

  ambientPeak <- round(median(cell_features_labeled[
    which(cell_features_labeled$barcode_class != "nucleus"),
  ]$num_transcripts, na.rm = TRUE))

  p5 <- ggdraw(function() {
    plotSelectedCellsSmoothScatter(cell_features_labeled,
      transcriptFeature = "num_transcripts",
      strTitlePrefix = paste0(
        "SVM nuclei method, ambient peak(UMIs) ",
        ambientPeak
      )
    )
  })

  p6 <- ggplot() +
    theme_void()
  if (useCellBenderFeatures) {
    p6 <- ggdraw(function() {
      plotSelectedCellsSmoothScatter(cell_features_labeled,
        transcriptFeature = "num_retained_transcripts",
        strTitlePrefix = "SVM nuclei method",
        useCellBenderFeatures = useCellBenderFeatures
      )
    })
  }

  return(list(
    cellbender = p1, initialization = p2, selected_nuclei = p3,
    nuclei_probability = p4, selected_nuclei_density = p5,
    selected_nuclei_density_rb = p6
  ))
}


################## INPUT VALIDATION

parseInputs <- function(
  cellFeaturesFile = NULL, dgeMatrixFile = NULL,
  optimusH5File = NULL
) {
  msg1 <- paste0(
    "Supply either the cellFeaturesFile and dgeMatrixFile, or ",
    "optimusH5File.  Supplying all 3 is ambiguous."
  )

  if (!is.null(cellFeaturesFile) & !is.null(dgeMatrixFile) &
    !is.null(optimusH5File)) {
    stop(msg1)
  }

  if (is.null(cellFeaturesFile) & is.null(dgeMatrixFile) &
    is.null(optimusH5File)) {
    msg <- paste0(
      "Supply either the cellFeaturesFile and dgeMatrixFile ",
      " or optimusH5File."
    )
    stop(msg)
  }

  if (!is.null(optimusH5File)) {
    r <- parseOptimusH5ad(optimusH5File)
    result <- list(cell_features = r$cell_features, dgeMatrix = r$dge)
    return(result)
  }


  cell_features <- readCellFeatures(cellFeaturesFile)
  dgeMatrix <- readDgeFile(dgeMatrixFile, cell_features)

  # filter out cell barcodes with expression values of 0.
  r <- filterZeroExpressionBarcodes(cell_features, dgeMatrix)
  cell_features <- r$cell_features
  dgeMatrix <- r$dgeMatrix

  result <- list(cell_features = cell_features, dgeMatrix = dgeMatrix)
  return(result)
}

filterZeroExpressionBarcodes <- function(cell_features, dgeMatrix) {
  # filter out cell barcodes with no expression
  idx <- which(Matrix::colSums(dgeMatrix) == 0)
  if (length(idx) > 0) {
    n <- min(5, length(idx))
    log_warn(
      "Removing cell barcodes with no expression ",
      "first 5 examples of [", length(idx), "] ",
      paste(head(colnames(dgeMatrix)[idx], n = n), collapse = ",")
    )
    dgeMatrix <- dgeMatrix[, -idx]
    cell_features <-
      cell_features[cell_features$cell_barcode %in% colnames(dgeMatrix), ]
  }

  return(list(cell_features = cell_features, dgeMatrix = dgeMatrix))
}

###############################################
# FIND EXEMPLAR CLASS BOUNDS
###############################################


#' Calculate the silhouette score for a clustering of cells.
#'
#' @param cell_features_labeled must include a column named
#'   'training_label_class' that contains the exemplar labels ('empty' and
#'   'nucleus'). NA entries and other labels are not included in the training data
#'   set.
#' @param downsampleRate The fraction of the data to use for silhouette
#'   calculation. Default is 0.1.
#' @param showPlot If TRUE, a silhouette plot will be displayed.
#' @param verbose If TRUE, the mean silhouette score will be printed to the log
#' @return The mean silhouette score for the clustering.
#' @importFrom cluster silhouette
#' @importFrom scales rescale
#' @noRd
calculate_silhouette <- function(
  cell_features_labeled, downsampleRate = 0.1,
  showPlot = FALSE, verbose = FALSE
) {
  if (!"training_label_class" %in% colnames(cell_features_labeled)) {
    if (verbose) {
      log_info("training_label_class not found. Cannot calculate silhouette.")
    }
    return(list(mean_silhouette = NA_real_, data = data.frame()))
  }

  idx <- which(cell_features_labeled$training_label_class %in%
    c("empty", "nucleus"))
  d <- cell_features_labeled[idx, ]

  if (nrow(d) < 2) {
    if (verbose) {
      log_info("Not enough training rows to calculate silhouette.")
    }
    return(list(mean_silhouette = NA_real_, data = d))
  }

  sample_size <- max(1, floor(nrow(d) * downsampleRate))
  d <- d[sample(nrow(d), sample_size), ]

  d$cluster_numeric <- as.numeric(factor(d$training_label_class))
  if (length(unique(d$cluster_numeric)) < 2) {
    if (verbose) {
      log_info("Not enough clusters to calculate silhouette.")
    }
    result <- list(mean_silhouette = NA_real_, data = d)
    return(result)
  }

  data_to_cluster <- d[, c("num_transcripts", "pct_intronic")]

  data_to_cluster <- data.frame(
    num_transcripts =
      scales::rescale(log10(data_to_cluster$num_transcripts + 1)),
    pct_intronic = scales::rescale(data_to_cluster$pct_intronic)
  )

  silhouette_scores <-
    cluster::silhouette(d$cluster_numeric, stats::dist(data_to_cluster))

  if (showPlot) {
    plot(silhouette_scores, main = "Silhouette plot")
  }

  data_to_cluster$silhouette_score <- silhouette_scores[, 3]
  data_to_cluster$cluster_numeric <- d$cluster_numeric

  mean_silhouette <- mean(data_to_cluster$silhouette_score)

  if (verbose) {
    log_info(sprintf("Mean silhouette score: %.3f", mean_silhouette))
  }

  result <- list(mean_silhouette = mean_silhouette, data = data_to_cluster)
  return(result)
}

#' Density bounds selection with enforced unimodal distribution via smoothing
#'
#' If the data is bimodal at the default bandwidth, the bandwidth is increased
#' until a unimodal distribution is found.
#'
#' @param cell_features A data frame containing the cell features.
#' @param yAxisFeature A string specifying the feature to use for the y-axis.
#' @param pctDensity The density percentile to use for the bounds.
#' @param maxPeaksExpected The maximum number of peaks expected in the data.
#'   The data will be iteratively smoothed until this many peaks detected.
#' @param showPlot A boolean indicating whether to show the plot - this may
#'   display multiple plots along the search space.
#' @import hdrcde
#' @return A data frame containing the bounds for the UMI and intronic features.
#' @noRd
getHighestDensityIntervalsEnforcedSmoothing <- function(
  cell_features,
  yAxisFeature = "pct_intronic", pctDensity = 50, maxPeaksExpected = 1,
  showPlot = FALSE
) {
  probList <- unique(c(25, 50, 75, 90, 95, 99, pctDensity))
  probList <- probList[probList <= pctDensity]

  umiBounds <- getBoundsByDensity(
    x = log10(cell_features$num_transcripts),
    probList = probList, pctDensity = pctDensity,
    maxPeaksExpected = maxPeaksExpected,
    showPlot = showPlot
  )$intervals
  yBounds <- getBoundsByDensity(
    x = cell_features[[yAxisFeature]],
    probList = probList, pctDensity = pctDensity,
    maxPeaksExpected = maxPeaksExpected, showPlot = showPlot
  )$intervals
  df <- data.frame(
    umi_lower_bound = umiBounds[1],
    umi_upper_bound = umiBounds[2], intronic_lower_bound = yBounds[1],
    intronic_upper_bound = yBounds[2]
  )
  return(df)
}


# Need the maximum density under each interval!  if expecting more
# than on peak, may need to break ties selecting the peak with the
# highest density.  it may be useful to override the bandwidth with
# the iterative bandwidth used for the entire experiment, instead of
# coming up with a new estimate here.
getBoundsByDensity <- function(
  x, probList, pctDensity, maxPeaksExpected = 1,
  showPlot = FALSE
) {
  pct <- paste0(pctDensity, "%")
  unimodal <- FALSE

  # the same bandwidth used by hdrcde::hdr.
  h <- hdrcde::hdrbw(hdrcde::BoxCox(x, 1), mean(probList))

  # multiply bandwidth by 2 each iteration
  mult <- 1
  while (!unimodal) {
    h <- h * mult
    zHDR <- hdrcde::hdr(x, prob = probList, h = h)$hdr
    if (showPlot) {
      hdrcde::hdr.den(x, prob = probList, h = h)
    }

    result <- zHDR[pct, ]
    result <- result[!is.na(result)]
    if (length(result) <= 2 * maxPeaksExpected) {
      resultDF <- list(intervals = result, bw = h)
      return(resultDF)
    } else {
      mult <- mult * 2
    }
  }
}


selectNucleiExemplarBounds <- function(
  cell_features,
  maxContaminationThreshold = 0.1, max_umis_empty = 50, initialDensity = 95,
  bounds_empty = NULL, extendCellSelectionBounds = TRUE
) {
  # the interval for cell barcodes that are not empty for training
  # data.  set a threshold of the lowest 25% of contamination for
  # cells that have contamination < 1.
  cell_features_non_empty <- cell_features[cell_features$frac_contamination <
    1, ]
  contaminationThreshold <- stats::quantile(
    cell_features_non_empty$frac_contamination,
    probs = seq(0, 1, 0.01)
  )["25%"]

  cell_features_non_empty <-
    cell_features_non_empty[cell_features_non_empty$frac_contamination <=
      contaminationThreshold, ]
  bounds_non_empty <-
    getHighestDensityIntervalsEnforcedSmoothing(cell_features_non_empty,
      pctDensity = 95
    )

  # if requested, try to extend the initial selection to a larger
  # region.
  bounds_non_empty_extended <- NULL
  if (extendCellSelectionBounds & !is.null(bounds_empty)) {
    idx <- which(
      log10(cell_features$num_transcripts) >
        (bounds_empty$umi_lower_bound + 1) &
        cell_features$pct_intronic >=
          bounds_non_empty$intronic_lower_bound &
        cell_features$pct_intronic <=
          bounds_non_empty$intronic_upper_bound &
        cell_features$frac_contamination <=
          maxContaminationThreshold
    )

    bounds_non_empty_extended <-
      getHighestDensityIntervalsEnforcedSmoothing(cell_features[idx, ],
        pctDensity = 95, showPlot = FALSE
      )

    bounds_non_empty_extended <- merge_bounds(
      bounds_non_empty,
      bounds_non_empty_extended
    )
  }

  return(list(
    bounds_non_empty = bounds_non_empty,
    bounds_non_empty_extended = bounds_non_empty_extended
  ))
}


################################### USE BOUNDS TO LABEL TRAINING DATA

#' Label the training data based on the given bounds
#'
#' @param cell_features A data frame containing the cell features.
#' @param bounds_empty A data frame containing the bounds for the empty cells.
#' @param bounds_non_empty A data frame containing the bounds for the non-empty
#'   cells.
#' @param maxContaminationThreshold The maximum contamination threshold for
#'   non-empty cells (only when using cellbender features.)
#' @param useCellBenderFeatures A boolean indicating whether to use CellBender
#'   features.
#' @param verbose A boolean indicating whether to print log messages.
#' @return A data frame containing the training data with a new column
#'   `training_label_class`.
#' @noRd
labelTrainingData <- function(
  cell_features, bounds_empty, bounds_non_empty,
  maxContaminationThreshold = 0.1, useCellBenderFeatures = TRUE,
  training_empty_barcodes = NULL, training_nucleus_barcodes = NULL,
  training_debris_barcodes = NULL, bounds_debris = NULL,
  verbose = TRUE
) {
  if (useCellBenderFeatures) {
    result <- (labelTrainingDataCBRB(cell_features, bounds_empty,
      bounds_non_empty,
      maxContaminationThreshold = maxContaminationThreshold,
      verbose = verbose
    ))
  } else {
    result <- labelTrainingDataDefault(cell_features, bounds_empty,
      bounds_non_empty,
      training_empty_barcodes = training_empty_barcodes,
      training_nucleus_barcodes = training_nucleus_barcodes,
      training_debris_barcodes = training_debris_barcodes,
      bounds_debris = bounds_debris,
      verbose = verbose
    )
  }
  return(result)
}

logTrainingDataSelection <- function(cell_features_labeled) {
  numEmpty <- 0
  numNonEmpty <- 0
  numDebris <- 0
  if ("training_label_class" %in% colnames(cell_features_labeled)) {
    numEmpty <- sum(
      cell_features_labeled$training_label_class == "empty",
      na.rm = TRUE
    )
    numNonEmpty <- sum(
      cell_features_labeled$training_label_class == "nucleus",
      na.rm = TRUE
    )
    numDebris <- sum(
      cell_features_labeled$training_label_class == "debris",
      na.rm = TRUE
    )
  }
  log_info("Number of empty exemplars: [", numEmpty, "]")
  log_info("Number of nuclei exemplars: [", numNonEmpty, "]")
  log_info("Number of debris exemplars: [", numDebris, "]")
}

labelTrainingDataCBRB <- function(
  cell_features, bounds_empty, bounds_non_empty,
  bounds_non_empty_extended = NULL, maxContaminationThreshold = 0.1,
  verbose = TRUE
) {
  # Define a function to find indices based on the given bounds
  find_indices <-
    function(cell_features, bounds, maxContaminationThreshold = 1,
             minContaminationThreshold = 0) {
      which(log10(cell_features$num_transcripts) >= bounds$umi_lower_bound &
        log10(cell_features$num_transcripts) <= bounds$umi_upper_bound &
        cell_features$pct_intronic >= bounds$intronic_lower_bound &
        cell_features$pct_intronic <= bounds$intronic_upper_bound &
        cell_features$frac_contamination <= maxContaminationThreshold &
        cell_features$frac_contamination >= minContaminationThreshold)
    }

  # Get the indices for empty and non-empty classes empty classes
  # have no filtering on max contamination - they are just empty -
  # all UMIs always subtracted from the ambient.
  idxEmpty <- find_indices(cell_features, bounds_empty, 1, 1)
  # non empty classes have a max contamination threshold
  if (!is.null(bounds_non_empty_extended)) {
    idxNonEmpty <- find_indices(
      cell_features, bounds_non_empty_extended,
      maxContaminationThreshold, 0
    )
  } else {
    idxNonEmpty <- find_indices(
      cell_features, bounds_non_empty,
      maxContaminationThreshold, 0
    )
  }
  all <- sort(union(idxEmpty, idxNonEmpty))

  # Assign classes based on the indices
  training_data <- cell_features
  training_data$training_label_class <- NA_character_
  training_data$training_label_class[idxNonEmpty] <- "nucleus"
  training_data$training_label_class[idxEmpty] <- "empty"

  if (verbose) {
    log_info("Training empty/nuclei cell barcodes selected [using CBRB]")
  }
  return(training_data)
}

labelTrainingDataDefault <- function(
  cell_features,
  bounds_empty,
  bounds_non_empty,
  training_empty_barcodes = NULL,
  training_nucleus_barcodes = NULL,
  training_debris_barcodes = NULL,
  bounds_debris = NULL,
  verbose = TRUE
) {
  if (!is.null(training_empty_barcodes) &&
    !is.null(training_nucleus_barcodes)) {
    return(labelTrainingDataDefaultByBarcode(
      cell_features = cell_features,
      training_empty_barcodes = training_empty_barcodes,
      training_nucleus_barcodes = training_nucleus_barcodes,
      training_debris_barcodes = training_debris_barcodes,
      verbose = verbose
    ))
  }

  # Define a function to find indices based on the given bounds
  find_indices <- function(cell_features, bounds) {
    which(log10(cell_features$num_transcripts) >= bounds$umi_lower_bound &
      log10(cell_features$num_transcripts) <= bounds$umi_upper_bound &
      cell_features$pct_intronic >= bounds$intronic_lower_bound &
      cell_features$pct_intronic <= bounds$intronic_upper_bound)
  }

  # Get the indices for empty and non-empty classes empty classes
  # have no filtering on max contamination - they are just empty -
  # all UMIs always subtracted from the ambient.
  idxEmpty <- find_indices(cell_features, bounds_empty)
  # non empty classes have a max contamination threshold
  idxNonEmpty <- find_indices(cell_features, bounds_non_empty)
  idxDebris <- integer(0)
  if (!is.null(bounds_debris) && !any(is.na(bounds_debris))) {
    idxDebris <- find_indices(cell_features, bounds_debris)
  }
  all <- sort(union(union(idxEmpty, idxNonEmpty), idxDebris))

  # Assign classes based on the indices
  training_data <- cell_features
  training_data$training_label_class <- NA_character_
  training_data$training_label_class[idxNonEmpty] <- "nucleus"
  training_data$training_label_class[idxEmpty] <- "empty"
  training_data$training_label_class[idxDebris] <- "debris"

  if (verbose) {
    log_info(
      "Training empty/nuclei cell barcodes selected",
      "[using density-only method]"
    )
  }
  return(training_data)
}

labelTrainingDataDefaultByBarcode <- function(
  cell_features,
  training_empty_barcodes,
  training_nucleus_barcodes,
  training_debris_barcodes = NULL,
  verbose = TRUE
) {
  if (is.null(rownames(cell_features)) || any(rownames(cell_features) == "")) {
    stop("Cell feature row names are required for barcode-based training labels.")
  }

  training_data <- cell_features
  training_data$training_label_class <- NA_character_

  idxEmpty <- which(rownames(training_data) %in% training_empty_barcodes)
  idxNonEmpty <- which(rownames(training_data) %in% training_nucleus_barcodes)
  idxDebris <- which(rownames(training_data) %in% training_debris_barcodes)

  training_data$training_label_class[idxNonEmpty] <- "nucleus"
  training_data$training_label_class[idxEmpty] <- "empty"
  training_data$training_label_class[idxDebris] <- "debris"

  if (verbose) {
    log_info(
      "Training empty/nuclei cell barcodes selected",
      "[using barcode-defined training labels]"
    )
  }

  training_data
}

# Function to merge dataframes
merge_bounds <- function(df1, df2) {
  merged_df <- df1
  for (col in names(df1)) {
    if (grepl("lower", col)) {
      merged_df[[col]] <- pmin(df1[[col]], df2[[col]])
    } else if (grepl("upper", col)) {
      merged_df[[col]] <- pmax(df1[[col]], df2[[col]])
    }
  }
  return(merged_df)
}

###################################### RUN SVM

runBinarySVM <- function(
  cell_features_labeled,
  features,
  positive_class,
  negative_class,
  probability_col,
  min_positive_exemplars = 100,
  min_negative_exemplars = 100
) {
  if (!"training_label_class" %in% colnames(cell_features_labeled)) {
    logger::log_warn(
      "training_label_class not found. Skipping SVM [",
      probability_col, "]."
    )
    return(NULL)
  }

  num_positive <- sum(
    cell_features_labeled$training_label_class == positive_class,
    na.rm = TRUE
  )
  num_negative <- sum(
    cell_features_labeled$training_label_class == negative_class,
    na.rm = TRUE
  )

  if (num_positive < min_positive_exemplars ||
    num_negative < min_negative_exemplars) {
    logger::log_warn(
      "Skipping SVM [", probability_col, "] because class counts are ",
      positive_class, " [", num_positive, "] and ",
      negative_class, " [", num_negative, "]"
    )
    return(NULL)
  }

  missing_features <- setdiff(features, colnames(cell_features_labeled))
  if (length(missing_features) > 0) {
    logger::log_warn(
      "Skipping SVM [", probability_col, "] because feature(s) were missing: ",
      paste(missing_features, collapse = ", ")
    )
    return(NULL)
  }

  unusable_features <- features[vapply(features, function(feature) {
    values <- cell_features_labeled[[feature]]
    is.numeric(values) && !any(is.finite(values))
  }, logical(1))]
  if (length(unusable_features) > 0) {
    logger::log_warn(
      "Skipping SVM [", probability_col, "] because feature(s) had no finite values: ",
      paste(unusable_features, collapse = ", ")
    )
    return(NULL)
  }

  cell_features_scaled <- scaleFeatures(cell_features_labeled, features)

  training_idx <- which(
    cell_features_scaled$training_label_class %in%
      c(positive_class, negative_class)
  )

  trainingData <- cell_features_scaled[
    training_idx,
    c(features, "training_label_class")
  ]

  trainingData$training_label_class <- factor(
    trainingData$training_label_class,
    levels = c(negative_class, positive_class)
  )

  svm_model <- e1071::svm(training_label_class ~ .,
    data = trainingData,
    type = "C-classification", kernel = "radial", probability = TRUE
  )

  predictions <- stats::predict(
    svm_model,
    cell_features_scaled[, features, drop = FALSE],
    probability = TRUE
  )

  probabilities <- attr(predictions, "probabilities")

  list(
    probabilities = probabilities[, positive_class],
    trainingData = trainingData,
    svm_model = svm_model
  )
}

scaleFeatures <- function(cell_features_labeled, features) {
  # make an explicity copy
  cell_features_scaled <- data.frame(cell_features_labeled)
  # log 10 the larger features
  cell_features_scaled$num_transcripts <-
    log10(cell_features_scaled$num_transcripts)

  # This feature isn't used by default but is scaled in case
  # someone wants to experiment.
  if ("num_reads" %in% cell_features_scaled) {
    cell_features_scaled$num_reads <- log10(cell_features_scaled$num_reads)
  }

  # scale the features prior to SVN
  cell_features_scaled[, features] <- scale(cell_features_scaled[, features])

  return(cell_features_scaled)
}


######################## PLOTTING CODE


#' Create a single page of cell selection plots
#'
#' @param plots A list of ggplot2 plots generated by callByIntronicSVM
#' @param geneModulePlots (Optional) A list of ggplot2 plots that capture the
#'   expression feautres of the gene modules
#' @param featurePlot (Optional) A ggplot2 plot that captures the feature values
#'   for the training data.
#' @param dataset_name The name of the dataset
#' @param outPDF The PDF file to write the plots to
#' @param useOpenPDF If TRUE, do not open/close the outPDF.  If false, opens a
#'   new PDF device, writes the plots, and closes the device.
#' @import cowplot ggplot2
#' @noRd
arrangeSVMCellSelectionPlots <- function(
  plots,
  geneModulePlots = NULL,
  featurePlotNucleusVsEmpty = NULL,
  featurePlotNucleusVsDebris = NULL,
  geneModuleExemplarPlot = NULL,
  dataset_name,
  outPDF,
  useOpenPDF = FALSE
) {
  # plots 1, 3, and 4 are ggplot2 objects.
  plots[[1]] <- plots[[1]] + custom_theme()
  plots[[3]] <- plots[[3]] + custom_theme()
  plots[[4]] <- plots[[4]] + custom_theme()

  plot_grid <- cowplot::plot_grid(
    plotlist = plots,
    ncol = 2,
    nrow = 3
  )

  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      dataset_name,
      fontface = "bold",
      size = 12,
      hjust = 0.5
    )

  final_plot <- cowplot::plot_grid(
    title,
    plot_grid,
    ncol = 1,
    rel_heights = c(0.05, 1)
  )

  if (!useOpenPDF && !is.null(outPDF)) {
    grDevices::pdf(outPDF)
  }

  # Page 1: main CBRB selection summary.
  gridExtra::grid.arrange(final_plot)

  # Page 2: selected-cell probability plot with larger text.
  plots[[4]] <- plots[[4]] +
    custom_theme(
      title_size = 12,
      axis_title_size = 10,
      axis_text_size = 10,
      legend_title_size = 10,
      legend_text_size = 10
    )

  gridExtra::grid.arrange(plots[[4]])

  # Page 3: nucleus vs empty feature plot.
  if (!is.null(featurePlotNucleusVsEmpty)) {
    gridExtra::grid.arrange(featurePlotNucleusVsEmpty)
  }

  # Page 4: nucleus vs debris feature plot, if the debris SVM ran.
  if (!is.null(featurePlotNucleusVsDebris)) {
    gridExtra::grid.arrange(featurePlotNucleusVsDebris)
  }

  # Page 5: debris gene module diagnostics, if the debris SVM ran.
  if (!is.null(geneModuleExemplarPlot) &&
    !is.null(geneModulePlots[["debris_gene_module_score"]])) {
    gridExtra::grid.arrange(
      geneModulePlots[["debris_gene_module_score"]] +
        custom_theme(),
      geneModuleExemplarPlot,
      ncol = 1
    )
  }

  # Existing gene-module summary page.
  if (!is.null(geneModulePlots)) {
    gridExtra::grid.arrange(
      arrangeSVMGeneModulePlots(
        geneModulePlots,
        dataset_name
      )
    )
  }

  if (!useOpenPDF && !is.null(outPDF)) {
    grDevices::dev.off()
  }
}

#' Create set of cell selection plots to evaluate cell selectionwhen CBRB
#' features are not used.
#'
#' @param plots A list of ggplot2 plots generated by callByIntronicSVM
#' @param geneModulePlots (Optional) A list of ggplot2 plots that capture the
#'   expression feautres of the gene modules
#' @param featurePlot (Optional) A ggplot2 plot that captures the feature values
#'   for the training data.
#' @param dataset_name The name of the dataset
#' @param outPDF The PDF file to write the plots to
#' @param useOpenPDF If TRUE, do not open/close the outPDF.  If false, opens a
#'   new PDF device, writes the plots, and closes the device.
#' @import cowplot ggplot2
#' @noRd
arrangeSVMCellSelectionPlotsNoCBRB <- function(
  plots,
  geneModulePlots = NULL,
  featurePlotNucleusVsEmpty = NULL,
  featurePlotNucleusVsDebris = NULL,
  geneModuleExemplarPlot = NULL,
  dataset_name,
  outPDF,
  useOpenPDF = FALSE
) {
  # re-order plots to include some plots that would otherwise be
  # empty due to CBRB features not being used.
  # 1. initial selection exemplars
  # 2. empty score distribution based on exemplars
  # 3. empty score distribution
  # 4. Selected cells
  # 5. Selected cell probability
  # 6. SmoothScatter final selection

  pList <- list(
    plots[[2]],
    geneModulePlots[["empty_gene_module_score_training_data"]] +
      custom_theme(),
    geneModulePlots[["empty_gene_module_score"]] +
      custom_theme(),
    plots[[3]] + custom_theme(),
    plots[[4]] + custom_theme(),
    plots[[5]]
  )

  plot_grid <- cowplot::plot_grid(
    plotlist = pList,
    ncol = 2,
    nrow = 3
  )

  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      dataset_name,
      fontface = "bold",
      size = 12,
      hjust = 0.5
    )

  final_plot <- cowplot::plot_grid(
    title,
    plot_grid,
    ncol = 1,
    rel_heights = c(0.05, 1)
  )

  if (!useOpenPDF && !is.null(outPDF)) {
    grDevices::pdf(outPDF)
  }

  gridExtra::grid.arrange(final_plot)

  if (!is.null(featurePlotNucleusVsEmpty)) {
    gridExtra::grid.arrange(featurePlotNucleusVsEmpty)
  }

  if (!is.null(featurePlotNucleusVsDebris)) {
    gridExtra::grid.arrange(featurePlotNucleusVsDebris)
  }

  if (!is.null(geneModuleExemplarPlot) &&
    !is.null(geneModulePlots[["debris_gene_module_score"]])) {
    gridExtra::grid.arrange(
      geneModulePlots[["debris_gene_module_score"]] +
        custom_theme(),
      geneModuleExemplarPlot,
      ncol = 1
    )
  }

  if (!useOpenPDF && !is.null(outPDF)) {
    grDevices::dev.off()
  }
}

#' Create a single page of gene module plots.
#'
#' @param plots A list of ggplot2 plots generated by callByIntronicSVM
#' @param dataset_name The name of the dataset
#' @import cowplot ggplot2
#' @noRd
arrangeSVMGeneModulePlots <- function(plots, dataset_name) {
  # all plots are ggplot2 plots.

  plots_custom_theme <- lapply(plots, function(plot) plot + custom_theme())
  # https://github.com/satijalab/seurat/issues/3584
  plots_custom_theme[["training_data"]] <- plots[["training_data"]] &
    custom_theme()

  # final plots to include:
  plotList <- c(
    "empty_gene_module_score_training_data",
    "empty_gene_module_score",
    "debris_gene_module_score_training_data",
    "debris_gene_module_score",
    "empty_gene_module_score_vs_contam",
    "cell_probability_histogram"
  )
  plotList <- plotList[plotList %in% names(plots_custom_theme)]
  plots_custom_theme <- plots_custom_theme[plotList]

  # Arrange the plots into a grid
  nrow <- ceiling(length(plots_custom_theme) / 2)
  plot_grid <- plot_grid(plotlist = plots_custom_theme, ncol = 2, nrow = nrow)

  # Add the title
  title <- ggdraw() + draw_label(dataset_name,
    fontface = "bold", size = 12,
    hjust = 0.5
  )

  # Combine the title and the plot grid
  final_plot <- plot_grid(title, plot_grid, ncol = 1, rel_heights = c(
    0.05,
    1
  ))
  return(final_plot)
}

# A custom theme for ggplot2 plots to make the text size appropriate
# when all plots are put together on a single page.
custom_theme <- function(
  title_size = 8, axis_title_size = 6,
  axis_text_size = 6, legend_title_size = 6, legend_text_size = 6
) {
  theme(
    plot.title = element_text(size = title_size, face = "bold"),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    legend.title = element_text(size = legend_title_size),
    legend.text = element_text(size = legend_text_size)
  )
}


#' Get the title for the cell selection plot
#'
#' @param cell_features_result The cell features data frame
#' @param strTitlePrefix The prefix for the title
#' @param transcriptFeature The feature to use for the number of transcripts
#' @return The title for the plot
#' @importFrom stats median
#' @importFrom scales comma_format
#' @noRd
getCellSelectionPlotTitle <- function(
  cell_features_result, strTitlePrefix = "",
  transcriptFeature = "num_transcripts"
) {
  selected <- cell_features_result[which(
    cell_features_result$barcode_class == "nucleus" &
      cell_features_result[[transcriptFeature]] > 0
  ), ]

  numSTAMPs <- dim(selected)[1]
  readsPerUMI <- NA

  minNumUmis <- min(selected[[transcriptFeature]])
  medianUMI <- round(median(selected[[transcriptFeature]]))
  minIntronic <- min(selected$pct_intronic)
  comma_formatter <- scales::label_comma()

  readsPerUMI <- mean(selected$num_reads / selected[[transcriptFeature]])

  strTitle <- sprintf(
    "%s, intronic>=%.2f\n%s Nuclei, %d+ UMIs, medUMIs %s",
    strTitlePrefix, minIntronic, comma_formatter(numSTAMPs), minNumUmis,
    comma_formatter(medianUMI)
  )

  # if the number of reads is available, add the average reads/UMI
  # to the title.
  if ("num_reads" %in% colnames(selected)) {
    strTitle <- sprintf(
      paste0(
        "%s, intronic >= %.2f\n%s Nuclei, %.1f reads/UMI,",
        " %d+ UMIs, medUMIs %s"
      ),
      strTitlePrefix, minIntronic, comma_formatter(numSTAMPs),
      readsPerUMI, minNumUmis, comma_formatter(medianUMI)
    )
  }

  return(strTitle)
}

#' Plot the expression vs intronic content of the cells, colored by the fraction
#' of UMIs removed by CellBender.
#'
#' @param cell_features A data frame containing the cell features.
#' @param densityCenters A list containing the density centers.
#' @param intronicThreshold The intronic threshold to plot.
#' @param title The title of the plot.
#' @param point_size The size of the points.
#' @param alpha The alpha value of the points.
#' @param useCellBenderFeatures If false, an empty plot is returned instead.
#' @return A ggplot object.
#' @import ggplot2
#' @noRd
plotExpressionVsIntronic <- function(
  cell_features, densityCenters = NULL,
  intronicThreshold = NULL, title = "All cell barcodes", point_size = 0.25,
  alpha = 0.25, useCellBenderFeatures = TRUE
) {
  if (!useCellBenderFeatures) {
    p <- ggplot() +
      theme_void() # empty plot
    return(p)
  }

  # TO MAKE R CMD CHECK HAPPY
  x <- y <- num_transcripts <- pct_intronic <- frac_contamination <- NULL

  p <- ggplot(cell_features, aes(
    x = log10(num_transcripts),
    y = pct_intronic, color = frac_contamination
  )) +
    ggrastr::rasterize(geom_point(size = point_size, alpha = alpha), dpi = 900) +
    scale_color_gradientn(
      colors = c("light blue", "blue", "red"), values = c(0, 0.5, 1)
    ) +
    labs(
      x = "log10(UMI)",
      y = "% Intronic",
      color = "CellBender\nUMIs\nRemoved"
    ) +
    coord_cartesian(xlim = log10_UMI_AXIS_RANGE_NEW) +
    ggtitle(title) +
    theme_minimal() +
    custom_theme(title_size = 8, axis_title_size = 8, axis_text_size = 8)

  if (!is.null(densityCenters)) {
    p <- p +
      geom_point(
        data = data.frame(
          x = log10(densityCenters$empty[1]),
          y = densityCenters$empty[2]
        ),
        aes(x = x, y = y), color = "red", size = 3
      ) +
      geom_point(
        data = data.frame(
          x = log10(densityCenters$cell[1]),
          y = densityCenters$cell[2]
        ),
        aes(x = x, y = y), color = "red", size = 3
      )
  }

  if (!is.null(intronicThreshold)) {
    p <- p +
      geom_hline(
        yintercept = intronicThreshold, linetype = "dotted",
        color = "red", linewidth = 1
      )
  }
  return(p)
}


#' Plot the training exemplar selection intervals
#'
#' This function plots log10 UMI count versus pct_intronic for all cell barcodes
#' and overlays the selected empty-droplet and nucleus exemplar regions. The
#' original one-dimensional HDR selections are shown as rectangular bounds when
#' `show_1d_bounds` is `TRUE`. If barcode lists from the two-dimensional HDR
#' refinement are supplied, their selected outer edges are drawn as convex hulls
#' when `show_2d_bounds` is `TRUE`.
#'
#' The two-dimensional hulls are used only for visualization. The actual
#' training labels are determined upstream from the selected barcode lists.
#'
#' @param cell_features A data frame containing barcode-level features. The
#'   data frame must contain `num_transcripts`, `pct_intronic`, and
#'   `training_label_class`. If two-dimensional bounds are plotted, row names
#'   must contain the cell barcode identifiers used in
#'   `training_empty_barcodes` and `training_nucleus_barcodes`.
#' @param bounds_empty A one-row data frame containing the rectangular
#'   empty-droplet exemplar bounds. Expected columns are `umi_lower_bound`,
#'   `umi_upper_bound`, `intronic_lower_bound`, and `intronic_upper_bound`.
#' @param bounds_non_empty A one-row data frame containing the rectangular
#'   nucleus exemplar bounds. Expected columns are `umi_lower_bound`,
#'   `umi_upper_bound`, `intronic_lower_bound`, and `intronic_upper_bound`.
#' @param training_empty_barcodes Optional character vector of empty-droplet
#'   exemplar barcodes selected by the two-dimensional HDR refinement.
#' @param training_nucleus_barcodes Optional character vector of nucleus
#'   exemplar barcodes selected by the two-dimensional HDR refinement.
#' @param show_1d_bounds Logical scalar. If `TRUE`, draw the rectangular
#'   one-dimensional HDR bounds for empty droplets and nuclei.
#' @param show_2d_bounds Logical scalar. If `TRUE`, draw convex hulls around
#'   the supplied two-dimensional refined empty-droplet and nucleus exemplar
#'   barcodes.
#' @param strTitleOverride Optional character scalar. If supplied, use this as
#'   the plot title instead of the default training-exemplar summary.
#' @param cex.axis Numeric scalar controlling axis tick-label size.
#' @param cex.lab Numeric scalar controlling axis label size.
#' @param cex.main Numeric scalar controlling title size.
#'
#' @return Invisibly returns `NULL`. The function is called for its base R
#'   plotting side effect.
#' @importFrom graphics par smoothScatter title rect lines
#' @importFrom grDevices chull
#' @noRd
plotCellTypeIntervals <- function(
  cell_features, bounds_empty, bounds_non_empty,
  bounds_debris = NULL,
  training_empty_barcodes = NULL,
  training_nucleus_barcodes = NULL,
  training_debris_barcodes = NULL,
  umi_threshold = NULL,
  debris_intronic_floor = NULL,
  show_1d_bounds = TRUE,
  show_2d_bounds = TRUE,
  strTitleOverride = NULL, cex.axis = 0.6, cex.lab = 0.7, cex.main = 0.65
) {
  par(mar = c(2, 2, 2, 1), mgp = c(0.8, 0.25, 0), tck = -0.02)

  drawBarcodeHull <- function(
    cell_features,
    selected_barcodes,
    border,
    lwd = 2,
    lty = "13"
  ) {
    if (is.null(selected_barcodes) || length(selected_barcodes) < 3) {
      return(invisible(NULL))
    }

    if (is.null(rownames(cell_features))) {
      return(invisible(NULL))
    }

    idx <- rownames(cell_features) %in% selected_barcodes

    if (sum(idx) < 3) {
      return(invisible(NULL))
    }

    x <- log10(cell_features$num_transcripts[idx])
    y <- cell_features$pct_intronic[idx]

    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]

    if (length(x) < 3) {
      return(invisible(NULL))
    }

    hull_idx <- grDevices::chull(x, y)
    hull_idx <- c(hull_idx, hull_idx[1])

    graphics::lines(
      x[hull_idx],
      y[hull_idx],
      col = border,
      lwd = lwd,
      lty = lty
    )

    invisible(NULL)
  }

  num_training_empty <- 0
  num_training_nuclei <- 0
  num_training_debris <- 0
  if ("training_label_class" %in% colnames(cell_features)) {
    num_training_empty <- sum(
      cell_features$training_label_class == "empty",
      na.rm = TRUE
    )
    num_training_nuclei <- sum(
      cell_features$training_label_class == "nucleus",
      na.rm = TRUE
    )
    num_training_debris <- sum(
      cell_features$training_label_class == "debris",
      na.rm = TRUE
    )
  }

  graphics::smoothScatter(
    log10(cell_features$num_transcripts),
    cell_features$pct_intronic,
    xlim = log10_UMI_AXIS_RANGE_NEW,
    ylim = c(0, 1),
    xlab = "log10( UMIs )",
    ylab = "intronic",
    main = "",
    axes = TRUE,
    cex.axis = cex.axis,
    cex.lab = cex.lab
  )

  strTitle <- paste0(
    "SVM Training -- empty [",
    num_training_empty,
    "] nuclei [",
    num_training_nuclei,
    "] debris [",
    num_training_debris,
    "]"
  )

  if (!is.null(strTitleOverride)) {
    strTitle <- strTitleOverride
  }

  graphics::title(
    main = strTitle,
    line = 1,
    col.main = "black",
    cex.main = cex.main
  )

  if (show_1d_bounds) {
    graphics::rect(
      bounds_empty$umi_lower_bound,
      bounds_empty$intronic_lower_bound,
      bounds_empty$umi_upper_bound,
      bounds_empty$intronic_upper_bound,
      border = "red",
      lwd = 2,
      lty = 1
    )

    graphics::rect(
      bounds_non_empty$umi_lower_bound,
      bounds_non_empty$intronic_lower_bound,
      bounds_non_empty$umi_upper_bound,
      bounds_non_empty$intronic_upper_bound,
      border = "green",
      lwd = 2,
      lty = 1
    )

    if (!is.null(bounds_debris) && !any(is.na(bounds_debris))) {
      graphics::rect(
        bounds_debris$umi_lower_bound,
        bounds_debris$intronic_lower_bound,
        bounds_debris$umi_upper_bound,
        bounds_debris$intronic_upper_bound,
        border = "orange",
        lwd = 2,
        lty = 1
      )
    }
  }

  if (show_2d_bounds) {
    drawBarcodeHull(
      cell_features = cell_features,
      selected_barcodes = training_empty_barcodes,
      border = "red",
      lwd = 4,
      lty = 1
    )

    drawBarcodeHull(
      cell_features = cell_features,
      selected_barcodes = training_nucleus_barcodes,
      border = "green",
      lwd = 4,
      lty = 1
    )
  }

  if (!show_1d_bounds && !is.null(bounds_debris) &&
    !any(is.na(bounds_debris))) {
    graphics::rect(
      bounds_debris$umi_lower_bound,
      bounds_debris$intronic_lower_bound,
      bounds_debris$umi_upper_bound,
      bounds_debris$intronic_upper_bound,
      border = "orange",
      lwd = 2,
      lty = 1
    )
  }


  if (!is.null(umi_threshold) && !is.na(umi_threshold)) {
    graphics::abline(v = umi_threshold, lty = 2, col = "red")
  }

  if (!is.null(debris_intronic_floor) && !is.na(debris_intronic_floor)) {
    graphics::abline(h = debris_intronic_floor, lty = 2, col = "blue")
  }

  par(mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0), tck = NA)
}


plotSelectedCells <- function(cell_features_result, size = 0.25, alpha = 0.25) {
  strTitle <- "Selected Nuclei"

  df <- cell_features_result[cell_features_result$barcode_class == "nucleus", ]
  umi_min_threshold <- min(log10(df[["num_transcripts"]]))
  intronic_min_threshold <- min(df[["pct_intronic"]])

  # TO MAKE R CMD CHECK HAPPY
  num_transcripts <- pct_intronic <- barcode_class <- NULL

  p <- ggplot(cell_features_result, aes(
    x = log10(num_transcripts),
    y = pct_intronic, color = barcode_class
  )) +
    ggrastr::rasterize(geom_point(size = size, alpha = alpha), dpi = 900) +
    labs(
      x = "log10(UMI)",
      y = "% Intronic",
      color = "Barcode class"
    ) +
    ggtitle(strTitle) +
    scale_color_manual(
      values = c(
        "nucleus" = "green",
        "debris" = "orange",
        "empty_or_other" = "lightblue"
      ),
      na.value = "grey80"
    ) +
    coord_cartesian(xlim = log10_UMI_AXIS_RANGE_NEW) +
    theme_minimal() +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

  return(p)
}

# changePar makes this plot compatable with the standard analysis
# cell selection pipeline that uses this plot in-line with standard
# graphics.

#' Plot Selected Cells Smooth Scatter
#'
#' This function creates a smooth scatter plot for selected cells, showing the
#' relationship between transcript counts and percent intronic features.
#'
#' @param cell_features A dataframe containing cell metadata with at least the
#'   columns specified by `transcriptFeature` and `pct_intronic`.
#' @param strTitlePrefix A prefix for the plot title.
#' @param transcriptFeature The column name for transcript counts (default:
#'   'num_transcripts').
#' @param changePar Logical; if TRUE, adjusts graphical parameters for the plot.
#' @param useCellBenderFeatures Logical; if FALSE, returns an empty ggplot
#'   object.
#' @return A plot object or NULL.
#' @import ggplot2
#' @importFrom graphics par smoothScatter title abline
#' @noRd
plotSelectedCellsSmoothScatter <- function(
  cell_features, strTitlePrefix = "",
  transcriptFeature = "num_transcripts", changePar = TRUE,
  useCellBenderFeatures = TRUE
) {
  if (!useCellBenderFeatures) {
    p <- ggplot2::ggplot() +
      ggplot2::theme_void() # empty plot
    return(p)
  }

  xlab <- "log10 ( UMIs )"

  if (transcriptFeature == "num_retained_transcripts") {
    xlab <- "log10 ( UMIs post CBRB )"
  }

  cell_features_filtered <- cell_features[!is.na(cell_features$barcode_class), ]
  # for cases where plotting remove background processed data.
  # number of transcripts should never be 0.
  cell_features_filtered <-
    cell_features_filtered[cell_features_filtered[[transcriptFeature]] > 0, ]

  strTitle <- getCellSelectionPlotTitle(cell_features_filtered,
    strTitlePrefix = strTitlePrefix, transcriptFeature = transcriptFeature
  )
  df <- cell_features_filtered[cell_features_filtered$barcode_class == "nucleus" &
    cell_features_filtered[[transcriptFeature]] > 0, ]
  umi_min_threshold <- min(df[[transcriptFeature]])
  intronic_min_threshold <- min(df[["pct_intronic"]])

  if (changePar) {
    opar <- graphics::par(no.readonly = TRUE)
    graphics::par(mar = c(2, 2, 2, 1), mgp = c(0.8, 0.25, 0), tck = -0.02)
  }

  graphics::smoothScatter(log10(cell_features[[transcriptFeature]]),
    cell_features$pct_intronic,
    xlim = log10_UMI_AXIS_RANGE_NEW, ylim = c(
      0,
      1
    ), xlab = xlab, ylab = "intronic", axes = TRUE, cex.axis = 0.6,
    cex.lab = 0.7
  )

  graphics::title(
    main = strTitle, line = 0.25,
    col.main = "black", cex.main = 0.65
  )
  graphics::abline(v = log10(umi_min_threshold), col = "red", lty = 2)
  graphics::abline(h = intronic_min_threshold, col = "red", lty = 2)

  if (changePar) {
    graphics::par(opar)
  }
}

plotScaledTrainingDataFeatures <- function(
  trainingData,
  strTitle = "Training Data Features"
) {
  if (is.null(trainingData) || nrow(trainingData) == 0) {
    return(NULL)
  }

  label_col_index <- which(names(trainingData) == "training_label_class")
  if (length(label_col_index) != 1) {
    return(NULL)
  }

  feature_columns <- trainingData[, -label_col_index, drop = FALSE]
  label_column <- trainingData[, label_col_index]

  longData <- data.frame(
    feature = rep(names(feature_columns), each = nrow(trainingData)),
    value = unlist(feature_columns),
    training_label_class = rep(label_column, times = ncol(feature_columns))
  )

  feature <- value <- training_label_class <- NULL

  p <- ggplot(longData, aes(
    x = feature,
    y = value,
    fill = training_label_class
  )) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    scale_fill_manual(values = c(
      "nucleus" = "green",
      "empty" = "red",
      "debris" = "orange"
    )) +
    labs(
      title = strTitle,
      x = "Features",
      y = "Values",
      fill = "Training class"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    ) +
    coord_flip()

  return(p)
}

combineTrainingFeaturePlots <- function(...) {
  plot_list <- list(...)
  plot_list <- plot_list[!vapply(plot_list, is.null, logical(1))]

  if (length(plot_list) == 0) {
    return(NULL)
  }
  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  }

  cowplot::plot_grid(plotlist = plot_list, ncol = 1)
}


plotCellProbabilities <- function(
  cell_features,
  strTitle = "Nuclei Probability"
) {
  df <- cell_features[which(cell_features$is_cell_prob >= 0.5), ]

  breaks <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1)
  colors <- c("red", "purple", "orange", "blue", "green")

  # Assign bins to is_cell_prob
  df$is_cell_prob_bin <-
    cut(df$is_cell_prob, breaks = breaks, include.lowest = TRUE)

  # Plot

  # TO MAKE R CMD CHECK HAPPY
  num_transcripts <- pct_intronic <- is_cell_prob_bin <- NULL

  p <- ggplot(df, aes(
    x = log10(num_transcripts), y = pct_intronic,
    color = is_cell_prob_bin
  )) +
    ggrastr::rasterize(geom_point(size = 0.5, alpha = 0.25), dpi = 900) +
    labs(
      x = "log10(UMI)",
      y = "% Intronic",
      color = "Cell Probability"
    ) +
    # coord_cartesian(xlim = log10_UMI_AXIS_RANGE_NEW) +
    ggtitle(strTitle) +
    theme_minimal() +
    scale_color_manual(values = colors) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
  return(p)
}


####################### GENE MODULE PLOTS

scatterPlotModuleScore <- function(
  cell_features,
  moduleName = "nuclei_gene_module_score", strTitle = "",
  point_size = 0.25, alpha = 0.25
) {
  # TO MAKE R CMD CHECK HAPPY
  num_transcripts <- pct_intronic <- NULL

  if (!moduleName %in% colnames(cell_features)) {
    stop("Module name ", moduleName, " not found in cell features")
  }

  p <- ggplot(cell_features, aes(
    x = log10(num_transcripts),
    y = pct_intronic, color = get(moduleName)
  )) +
    ggrastr::rasterize(geom_point(
      size = point_size, alpha = alpha
    ), dpi = 900) +
    scale_color_gradient2(
      low = "red", mid = "lightgrey", high = "blue", midpoint = 0
    ) +
    labs(
      x = "log10(UMI)",
      y = "% Intronic",
      color = ""
    ) +
    coord_cartesian(xlim = log10_UMI_AXIS_RANGE_NEW) +
    ggtitle(moduleName) +
    theme_minimal() # +
  return(p)
}

plotModuleVsFracContam <- function(
  cell_features,
  moduleName = "nuclei_gene_module_score", point_size = 0.5, alpha = 0.5
) {
  # TO MAKE R CMD CHECK HAPPY
  frac_contamination <- pct_intronic <- NULL

  if (!moduleName %in% colnames(cell_features)) {
    stop("Module name ", moduleName, " not found in cell features")
  }

  p <- ggplot(cell_features, aes(
    x = frac_contamination, y = get(moduleName),
    color = pct_intronic
  )) +
    ggrastr::rasterize(geom_point(size = point_size, alpha = alpha),
      dpi = 900
    ) +
    labs(
      x = "Fraction Contamination",
      y = moduleName,
      color = "% Intronic"
    ) +
    theme_minimal()
  return(p)
}

plotCellProbabilityConditionalCbrb <- function(
  cell_features_result,
  frac_contamination_threshold = 1,
  useCellBenderFeatures = TRUE
) {
  # if not using celbender, short circuit with an empty plot
  if (!useCellBenderFeatures) {
    return(ggplot() +
      theme_minimal())
  }

  df <- cell_features_result[
    cell_features_result$frac_contamination >= frac_contamination_threshold,
  ]

  # Filter out NAs
  df <- df[!is.na(df$is_cell_prob), ]

  # Calculate the number of cell barcodes
  numCellBarcodes <- dim(df[df$barcode_class == "nucleus", ])[1]

  r <- graphics::hist(df$is_cell_prob,
    breaks = seq(0, 1, by = 0.01),
    plot = FALSE
  )

  dd <- data.frame(countLog10 = log10(r$counts + 1), midpoint = r$mids)

  # TO MAKE R CMD CHECK HAPPY
  midpoint <- countLog10 <- NULL

  title <- paste(
    "SVM Probability for CBRB empty cell barcodes\n[",
    numCellBarcodes, "] barcode selected as nuclei"
  )

  p <- ggplot(dd, aes(x = midpoint, y = countLog10)) +
    geom_bar(stat = "identity", width = 0.01) +
    geom_vline(
      xintercept = 0.5, color = "red",
      linetype = "dashed", linewidth = 1
    ) +
    labs(x = "Cell Probability", y = "Log10(Count + 1)") +
    ggtitle(title) +
    theme_minimal() +
    annotation_logticks(sides = "l")
  return(p)
}


#' Plot gene module scores for training exemplar classes
#'
#' This diagnostic plot compares the empty-droplet gene module score against the
#' debris gene module score for barcodes used as training exemplars. Points are
#' colored by `training_label_class`, so the plot shows whether the empty,
#' nucleus, and debris exemplar classes separate in the two-dimensional
#' gene-module score space.
#'
#' @param cell_features_labeled A cell-feature data frame containing
#'   `empty_gene_module_score`, `debris_gene_module_score`, and
#'   `training_label_class`.
#' @param strTitle Character scalar used as the plot title.
#' @param point_alpha Numeric alpha value used for points. Smaller values reduce
#'   overplotting.
#' @param point_size Numeric point size passed to [ggplot2::geom_point()].
#'
#' @return A ggplot object.
#'
#' @export
plotGeneModuleScoresByExemplarClass <- function(
  cell_features_labeled,
  strTitle = "Gene module scores by exemplar class",
  point_alpha = 0.25,
  point_size = 0.6
) {
  requiredCols <- c(
    "empty_gene_module_score",
    "debris_gene_module_score",
    "training_label_class"
  )

  missingCols <- setdiff(requiredCols, colnames(cell_features_labeled))
  if (length(missingCols) > 0) {
    stop(
      "Missing required column(s): ",
      paste(missingCols, collapse = ", ")
    )
  }

  idxKeep <- is.finite(cell_features_labeled$empty_gene_module_score) &
    is.finite(cell_features_labeled$debris_gene_module_score) &
    !is.na(cell_features_labeled$training_label_class)

  plotDf <- cell_features_labeled[idxKeep, , drop = FALSE]
  plotDf$training_label_class <- factor(
    plotDf$training_label_class,
    levels = c("empty", "nucleus", "debris")
  )

  # Make R CMD CHECK Happy
  empty_gene_module_score <- debris_gene_module_score <-
    training_label_class <- NULL

  ggplot2::ggplot(
    plotDf,
    ggplot2::aes(
      x = empty_gene_module_score,
      y = debris_gene_module_score,
      color = training_label_class
    )
  ) +
    ggplot2::geom_point(
      alpha = point_alpha,
      size = point_size
    ) +
    ggplot2::scale_color_manual(
      values = c(
        empty = "red",
        nucleus = "green",
        debris = "orange"
      ),
      drop = FALSE
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(
          size = point_size * 4,
          alpha = 1
        )
      )
    ) +
    ggplot2::labs(
      title = strTitle,
      x = "Empty gene module score",
      y = "Debris gene module score",
      color = "Exemplar class"
    ) +
    ggplot2::theme_classic()
}
