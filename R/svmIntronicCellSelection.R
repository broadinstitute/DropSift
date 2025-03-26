# Some global parameters
MIN_UMIs_PER_STAMP <- 20  # what's worth even including in analysis
log10_UMI_AXIS_RANGE_NEW <- c(log10(MIN_UMIs_PER_STAMP), 6)  # for plotting
random.seed <- 1

#' Call nuclei using an SVM trained on cell summary features.
#'
#' A convenience method for calling cells using an SVM trained on cell summary
#' features.
#'
#' Please set a random seed for reproducibility, e.g., `set.seed(1)`. This
#' ensures that the SVM initialization and training are consistent across runs.
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
#' @param forceTwoClusterSolution When true, the initialization of the SVM will
#'   attempt to find a solution with two clusters. In cases where an experiment
#'   is overloaded with nuclei, this may correct the initial set of nuclei and
#'   empty droplets selected. This argument is specific to the density-style
#'   selection of exemplars that is only applicable when useCBRBFeatures is
#'   false.
#' @param outPDF The PDF file to write the plots to.
#' @param outFeaturesFile The cell features dataframe, further annotated by the
#'   SVM to include the cell probability and label, along with which cell
#'   barcodes were used for training.
#' @param outCellBenderInitialParameters If useCBRBFeatures is false and this
#'   parameter is not null, the SVM is trained without cellbenber features, and
#'   the parameters expected_cells and total_droplets_included are
#' estimated and written to this file.
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
#'    datasetName = "example_dataset",
#'    cellFeaturesFile = cell_features_file,
#'    dgeMatrixFile = temp_dir,  # directory with compressed 10x files
#'    useCBRBFeatures = FALSE,
#'    forceTwoClusterSolution = FALSE,
#'    outPDF = out_pdf_file,
#'    outFeaturesFile = out_features_file
#')
#'
#' # Check that output files were generated
#' file.exists(out_pdf_file)    # Should return TRUE
#' file.exists(out_features_file)  # Should return TRUE
#'
#' # Inspect the PDF output report
#' if (file.exists(out_pdf_file)) {
#'    utils::browseURL(out_pdf_file)
#' }
runIntronicSVM <- function(datasetName, cellFeaturesFile = NULL,
    dgeMatrixFile = NULL, optimusH5File = NULL, cellProbabilityThreshold = NULL,
    max_umis_empty = 50, features = NULL, useCBRBFeatures = TRUE,
    forceTwoClusterSolution = FALSE, outPDF = NULL, outFeaturesFile = NULL,
    outCellBenderInitialParameters = NULL) {

    # parse and validate inputs.
    r <- parseInputs(cellFeaturesFile = cellFeaturesFile,
        dgeMatrixFile = dgeMatrixFile,
        optimusH5File = optimusH5File)

    cell_features <- r$cell_features
    dgeMatrix <- r$dgeMatrix
    svmNucleusCaller <- SvmNucleusCaller(cellFeatures = cell_features,
        dgeMatrix = dgeMatrix,
        cellProbabilityThreshold = cellProbabilityThreshold,
        maxUmisEmpty = max_umis_empty, featureColumns = features,
        forceTwoClusterSolution = forceTwoClusterSolution,
        useCBRBFeatures = useCBRBFeatures, datasetName = datasetName)

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
        write.table(cell_features_result, outFeaturesFile, row.names = FALSE,
            col.names = TRUE, quote = FALSE, sep = "\t")
    }

    # write initial cellbender remove-background parameters to file
    if (!is.null(outCellBenderInitialParameters)) {
        cbrbArgs <- getCBRBArgs(svmNucleusCaller)
        result <- data.frame(cbrbArgs)
        write.table(result, outCellBenderInitialParameters, row.names = FALSE,
            col.names = TRUE, quote = FALSE, sep = "\t")
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
#' @import e1071 hdrcde ggplot2 grid gridExtra cowplot ggrastr logger
#'   gridGraphics
#' @noRd
callByIntronicSVM <- function(dataset_name, cell_features, dgeMatrix,
    cellProbabilityThreshold = NULL, max_umis_empty = 50,
    features, forceTwoClusterSolution = FALSE) {

    validateFeaturePresence(cell_features, features)
    maxContaminationThreshold <- 0.1  # CBRB-specific contamination threshold
    log_info("Beginning SVM Intronic Cell Selection")
    featureListStr<-paste(features, collapse = ", ")
    log_info("Features used for cell selection:", featureListStr)
    rownames(cell_features) <- cell_features$cell_barcode
    cell_features$cell_barcode <- NULL
    useCellBenderFeatures <- "frac_contamination" %in% features

    allBounds <- findTrainingDataBounds(cell_features, max_umis_empty,
        useCellBenderFeatures = useCellBenderFeatures,
        forceTwoClusterSolution = forceTwoClusterSolution)
    bounds_empty <- allBounds$bounds_empty
    bounds_non_empty <- allBounds$bounds_non_empty
    cell_features_labeled <- labelTrainingData(cell_features, bounds_empty,
        bounds_non_empty, maxContaminationThreshold, useCellBenderFeatures)
    logTrainingDataSelection(cell_features_labeled)
    r <- addGeneModules(cell_features_labeled, dgeMatrix, numGenes = 100,
        useCellBenderFeatures = useCellBenderFeatures, verbose = FALSE)
    cell_features_labeled <- r$cell_features

    geneModulePlots <-
        r$gene_module_plots[seq_len(min(4, length(r$gene_module_plots)))]
    if (!"empty_gene_module_score" %in% colnames(cell_features_labeled)) {
        log_warn("empty_gene_module_score was not computed. Dropping feature.")
        features <- setdiff(features, "empty_gene_module_score")
    }

    svm_result <- runSVM(cell_features_labeled, features, bounds_empty,
        cellProbabilityThreshold)
    cell_features_result <- svm_result$cell_features
    trainingData <- svm_result$trainingData

    selectionPlots <- createSelectionVisualization(cell_features_result,
        bounds_empty, bounds_non_empty, useCellBenderFeatures, dataset_name)
    feature_plot <- plotScaledTrainingDataFeatures(
        svm_result$trainingData[, c(features, "training_label_is_cell")]
    )
    log_info("Nuclei selection finished")
    return(list(dataset_name = dataset_name,
        cell_features = svm_result$cell_features,
        features = features, plots = selectionPlots,
        featurePlot = feature_plot, geneModulePlots = geneModulePlots,
        bounds_empty = bounds_empty, bounds_non_empty = bounds_non_empty))
}

validateFeaturePresence <- function(cell_features, features) {
    for (f in setdiff(features, "empty_gene_module_score")) {
        if (!(f %in% colnames(cell_features))) {
            stop(paste("Requested Features [", paste(features, collapse = ", "),
                "]. Feature", f, "not found in cell_features data frame."))
        }
    }
}


createSelectionVisualization <- function(cell_features_labeled, bounds_empty,
    bounds_non_empty, useCellBenderFeatures, dataset_name) {
    p1 <- plotExpressionVsIntronic(cell_features_labeled,
        title = "All cell barcodes",
        useCellBenderFeatures = useCellBenderFeatures)

    p2 <- cowplot::ggdraw(function()
        plotCellTypeIntervals(cell_features_labeled, bounds_empty,
        bounds_non_empty))

    p3 <- plotSelectedCells(cell_features_labeled)
    p4 <- plotCellProbabilities(cell_features_labeled,
        strTitle = "Cell Probability")

    ambientPeak <- round(median(cell_features_labeled[
        which(cell_features_labeled$is_cell == FALSE), ]$num_transcripts))

    p5 <- ggdraw(function() {
        plotSelectedCellsSmoothScatter(cell_features_labeled,
            transcriptFeature = "num_transcripts",
            strTitlePrefix = paste0("SVM nuclei method, ambient peak(UMIs) ",
                ambientPeak))
    })

    p6 <- ggplot() + theme_void()
    if (useCellBenderFeatures) {
        p6 <- ggdraw(function() {
            plotSelectedCellsSmoothScatter(cell_features_labeled,
                transcriptFeature = "num_retained_transcripts",
                strTitlePrefix = "SVM nuclei method",
                useCellBenderFeatures = useCellBenderFeatures)
        })
    }

    return(list(cellbender = p1, initialization = p2, selected_nuclei = p3,
                nuclei_probability = p4, selected_nuclei_density = p5,
                selected_nuclei_density_rb = p6))
}


################## INPUT VALIDATION

parseInputs <- function(cellFeaturesFile = NULL, dgeMatrixFile = NULL,
    optimusH5File = NULL) {

    msg1 <- paste0("Supply either the cellFeaturesFile and dgeMatrixFile, or ",
        "optimusH5File.  Supplying all 3 is ambiguous.")

    if (!is.null(cellFeaturesFile) & !is.null(dgeMatrixFile)
        & !is.null(optimusH5File)) {
        stop(msg1)
    }

    if (is.null(cellFeaturesFile) & is.null(dgeMatrixFile)
        & is.null(optimusH5File)) {
        msg<- paste0("Supply either the cellFeaturesFile and dgeMatrixFile ",
            " or optimusH5File.")
        stop(msg)
    }

    if (!is.null(optimusH5File)) {
        r <- parseOptimusH5ad(optimusH5File)
        result <- list(cell_features = r$cell_features, dgeMatrix = r$dge)
        return(result)
    }


    cell_features <- readCellFeatures(cellFeaturesFile)
    dgeMatrix <- readDgeFile(dgeMatrixFile, cell_features)
    result <- list(cell_features = cell_features, dgeMatrix = dgeMatrix)
    return(result)
}

###############################################
#FIND EXEMPLAR CLASS BOUNDS
###############################################




#' Calculate the silhouette score for a clustering of cells.
#'
#' @param cell_features_labeled must include a column named
#'   'training_label_is_cell' that contains the cell/empty labels (TRUE for
#'   cell, FALSE for empty). NA entries are not included in the training data
#'   set.
#' @param downsampleRate The fraction of the data to use for silhouette
#'   calculation. Default is 0.1.
#' @param showPlot If TRUE, a silhouette plot will be displayed.
#' @param verbose If TRUE, the mean silhouette score will be printed to the log
#' @param seed The random seed to use for downsampling the data.
#' @return The mean silhouette score for the clustering.
#' @importFrom cluster silhouette
#' @importFrom scales rescale
#' @noRd
calculate_silhouette <- function(cell_features_labeled, downsampleRate = 0.1,
    showPlot = FALSE, verbose = FALSE) {
    # restrict to the training data
    idx <- which(!is.na(cell_features_labeled$training_label_is_cell))
    d <- cell_features_labeled[idx,]

    # downsample to 10% of the original data.
    d <- d[sample(nrow(d), nrow(d) * downsampleRate), ]

    # Ensure there are two clusters (TRUE and FALSE values)
    d$cluster_numeric <- as.numeric(factor(d$training_label_is_cell))
    if (length(unique(d$cluster_numeric)) < 2) {
        if (verbose)
            log_info("Not enough clusters to calculate silhouette.")
        result <- list(mean_silhouette = NA, data = d)
        return(result)
    }

    # Extract the features used for clustering
    data_to_cluster <- d[, c("num_transcripts", "pct_intronic")]

    # rescale the feaures to 0-1 scale
    data_to_cluster <- data.frame(
        num_transcripts <-
            scales::rescale(log10(data_to_cluster$num_transcripts +1)),
        pct_intronic <- scales::rescale(data_to_cluster$pct_intronic))

    # Compute the silhouette score based on the
    # training_label_is_cell column
    silhouette_scores <-
        cluster::silhouette(d$cluster_numeric, stats::dist(data_to_cluster))

    # Optionally show the silhouette plot
    if (showPlot) {
        plot(silhouette_scores, main = "Silhouette plot")
    }

    data_to_cluster$silhouette_score <- silhouette_scores[, 3]
    data_to_cluster$cluster_numeric <- d$cluster_numeric

    # Calculate the mean silhouette score for the clustering
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
getHighestDensityIntervalsEnforcedSmoothing <- function(cell_features,
    yAxisFeature = "pct_intronic", pctDensity = 50, maxPeaksExpected = 1,
    showPlot = FALSE) {
    probList <- unique(c(25, 50, 75, 90, 95, 99, pctDensity))
    probList <- probList[probList <= pctDensity]

    umiBounds <- getBoundsByDensity(x = log10(cell_features$num_transcripts),
        probList = probList, pctDensity = pctDensity,
        maxPeaksExpected = maxPeaksExpected,
        showPlot = showPlot)$intervals
    yBounds <- getBoundsByDensity(x = cell_features[[yAxisFeature]],
        probList = probList, pctDensity = pctDensity,
        maxPeaksExpected = maxPeaksExpected, showPlot = showPlot)$intervals
    df <- data.frame(umi_lower_bound = umiBounds[1],
        umi_upper_bound = umiBounds[2], intronic_lower_bound = yBounds[1],
        intronic_upper_bound = yBounds[2])
    return(df)
}


# Need the maximum density under each interval!  if expecting more
# than on peak, may need to break ties selecting the peak with the
# highest density.  it may be useful to override the bandwidth with
# the iterative bandwidth used for the entire experiment, instead of
# coming up with a new estimate here.
getBoundsByDensity <- function(x, probList, pctDensity, maxPeaksExpected = 1,
    showPlot = FALSE) {
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



selectNucleiExemplarBounds <- function(cell_features,
    maxContaminationThreshold = 0.1, max_umis_empty = 50, initialDensity = 95,
    bounds_empty = NULL, extendCellSelectionBounds = TRUE) {
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
        pctDensity = 95)

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
            getHighestDensityIntervalsEnforcedSmoothing(cell_features[idx,],
            pctDensity = 95, showPlot = FALSE)

        bounds_non_empty_extended <- merge_bounds(bounds_non_empty,
            bounds_non_empty_extended)
    }

    return(list(bounds_non_empty = bounds_non_empty,
                bounds_non_empty_extended = bounds_non_empty_extended))
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
#'   `training_label_is_cell`.
#' @noRd
labelTrainingData <- function(cell_features, bounds_empty, bounds_non_empty,
    maxContaminationThreshold = 0.1, useCellBenderFeatures = TRUE,
    verbose = TRUE) {

    if (useCellBenderFeatures) {
        result <- (labelTrainingDataCBRB(cell_features, bounds_empty,
            bounds_non_empty,
            maxContaminationThreshold = maxContaminationThreshold,
            verbose = verbose))
    } else {
        result <- labelTrainingDataDefault(cell_features, bounds_empty,
            bounds_non_empty, verbose = verbose)
    }
    return(result)
}

logTrainingDataSelection <- function(cell_features_labeled) {
    numEmpty <- sum(!is.na(cell_features_labeled$training_label_is_cell) &
        !cell_features_labeled$training_label_is_cell)
    numNonEmpty <- sum(!is.na(cell_features_labeled$training_label_is_cell) &
        cell_features_labeled$training_label_is_cell)
    log_info("Number of empty exemplars: [", numEmpty, "]")
    log_info("Number of nuclei exemplars: [", numNonEmpty, "]")
}

labelTrainingDataCBRB <- function(cell_features, bounds_empty, bounds_non_empty,
    bounds_non_empty_extended = NULL, maxContaminationThreshold = 0.1,
    verbose = TRUE) {
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
        idxNonEmpty <- find_indices(cell_features, bounds_non_empty_extended,
            maxContaminationThreshold, 0)
    } else {
        idxNonEmpty <- find_indices(cell_features, bounds_non_empty,
            maxContaminationThreshold, 0)
    }
    all <- sort(union(idxEmpty, idxNonEmpty))

    # Assign classes based on the indices
    training_data <- cell_features
    training_data$training_label_is_cell <- NA
    training_data$training_label_is_cell[idxNonEmpty] <- TRUE
    training_data$training_label_is_cell[idxEmpty] <- FALSE

    if (verbose) {
        log_info("Training empty/nuclei cell barcodes selected [using CBRB]")
    }
    return(training_data)
}

labelTrainingDataDefault <- function(cell_features, bounds_empty,
    bounds_non_empty, verbose = TRUE) {
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
    all <- sort(union(idxEmpty, idxNonEmpty))

    # Assign classes based on the indices
    training_data <- cell_features
    training_data$training_label_is_cell <- NA
    training_data$training_label_is_cell[idxNonEmpty] <- TRUE
    training_data$training_label_is_cell[idxEmpty] <- FALSE

    if (verbose) {
        log_info("Training empty/nuclei cell barcodes selected",
        "[using density-only method]")
    }
    return(training_data)
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

runSVM <- function(cell_features_labeled, features, bounds_empty,
    cellProbabilityThreshold = NULL) {
    # scale the requested cell features.
    cell_features_scaled <- scaleFeatures(cell_features_labeled, features)

    # get the training data - restrict data to the labeled training
    # data for the specific features we are using.
    trainingData <- cell_features_scaled[
        !is.na(cell_features_scaled$training_label_is_cell),
        c(features, "training_label_is_cell")
    ]

    trainingData$training_label_is_cell <-
        as.factor(trainingData$training_label_is_cell)

    # train the SVM
    svm_model <- e1071::svm(training_label_is_cell ~ ., data = trainingData,
        type = "C-classification", kernel = "radial", probability = TRUE)
    # predict all points
    predictions <- stats::predict(svm_model, cell_features_scaled[, features],
        probability = TRUE)
    # Extract the probabilities
    probabilities <- attr(predictions, "probabilities")
    # Combine predictions and confidence into a dataframe
    results_df <- data.frame(is_cell = predictions,
        is_cell_prob = probabilities[,"TRUE"])
    # optionally override cell probability
    if (!is.null(cellProbabilityThreshold)) {
        results_df$is_cell <-
            ifelse(results_df$is_cell_prob >= cellProbabilityThreshold,
            TRUE, FALSE)
    }

    # Exclude barcodes lower than the lower bound of UMIs from classification.
    idx <- which(log10(cell_features_labeled$num_transcripts) <
        bounds_empty$umi_upper_bound)

    results_df$is_cell_prob[idx] <- NA
    results_df$is_cell[idx] <- FALSE

    # merge the results_df with the original cell_features
    cell_features_result <- cbind(cell_features_labeled, results_df)
    result <- list(cell_features = cell_features_result,
        trainingData = trainingData, svm_model = svm_model)
    return(result)
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
arrangeSVMCellSelectionPlots <- function(plots, geneModulePlots = NULL,
    featurePlot = NULL, dataset_name, outPDF, useOpenPDF = FALSE) {
    # plots 1,3,4 are ggplot2.

    plots[[1]] <- plots[[1]] + custom_theme()
    plots[[3]] <- plots[[3]] + custom_theme()
    plots[[4]] <- plots[[4]] + custom_theme()

    # Arrange the plots into a grid
    plot_grid <- plot_grid(plotlist = plots, ncol = 2, nrow = 3)

    # Add the title
    title <- ggdraw() + draw_label(dataset_name, fontface = "bold", size = 12,
        hjust = 0.5)

    # Combine the title and the plot grid
    final_plot <- plot_grid(title, plot_grid, ncol = 1, rel_heights = c(0.05,
        1))

    # Set PDF dimensions pdf_width <- 12 pdf_height <- 12

    # Open the PDF device
    if (!useOpenPDF & !is.null(outPDF)) {
        # pdf(outPDF, width = pdf_width, height = pdf_height)
        grDevices::pdf(outPDF)
    }

    gridExtra::grid.arrange(final_plot)

    plots[[4]] <- plots[[4]] +
        custom_theme(title_size = 12, axis_title_size = 10,
        axis_text_size = 10, legend_title_size = 10, legend_text_size = 10)

    gridExtra::grid.arrange(plots[[4]])

    if (!is.null(featurePlot)) {
        gridExtra::grid.arrange(featurePlot)
    }

    if (!is.null(geneModulePlots)) {
        gridExtra::grid.arrange(arrangeSVMGeneModulePlots(geneModulePlots,
            dataset_name))
    }

    if (!useOpenPDF & !is.null(outPDF)) {
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
arrangeSVMCellSelectionPlotsNoCBRB <- function(plots, geneModulePlots = NULL,
    featurePlot = NULL, dataset_name, outPDF, useOpenPDF = FALSE) {
    # re-order plots to include some plots that would otherwise be
    # empty due to CBRB features not being used.
    # 1. initial selection exemplars
    # 2. empty score distribution based on exemplars
    # 3. empty score distribution
    # 4. Selected cells
    # 5. Selected cell probability
    # 6. SmoothScatter final selection

    pList <- list(
        plots[[2]], geneModulePlots[["training_data"]] + custom_theme(),
        geneModulePlots[["empty_gene_module_score"]] + custom_theme(),
        plots[[3]] + custom_theme(), plots[[4]] + custom_theme(), plots[[5]]
    )

    # Arrange the plots into a grid
    plot_grid <- plot_grid(plotlist = pList, ncol = 2, nrow = 3)

    # Add the title
    title <- ggdraw() + draw_label(dataset_name, fontface = "bold", size = 12,
        hjust = 0.5)

    # Combine the title and the plot grid
    final_plot <- plot_grid(title, plot_grid, ncol = 1, rel_heights = c(0.05,
        1))

    # Open the PDF device
    if (!useOpenPDF & !is.null(outPDF)) {
        # pdf(outPDF, width = pdf_width, height = pdf_height)
        grDevices::pdf(outPDF)
    }

    gridExtra::grid.arrange(final_plot)

    if (!is.null(featurePlot)) {
        gridExtra::grid.arrange(featurePlot)
    }

    if (!useOpenPDF & !is.null(outPDF)) {
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
    plotList <- c("training_data", "empty_gene_module_score",
        "empty_gene_module_score_vs_contam", "cell_probability_histogram")
    plots_custom_theme <- plots_custom_theme[plotList]

    # Arrange the plots into a grid
    nrow <- ceiling(length(plots_custom_theme)/2)
    plot_grid <- plot_grid(plotlist = plots_custom_theme, ncol = 2, nrow = nrow)

    # Add the title
    title <- ggdraw() + draw_label(dataset_name, fontface = "bold", size = 12,
        hjust = 0.5)

    # Combine the title and the plot grid
    final_plot <- plot_grid(title, plot_grid, ncol = 1, rel_heights = c(0.05,
        1))
    return(final_plot)
}

# A custom theme for ggplot2 plots to make the text size appropriate
# when all plots are put together on a single page.
custom_theme <- function(title_size = 8, axis_title_size = 6,
    axis_text_size = 6, legend_title_size = 6, legend_text_size = 6) {
    theme(plot.title = element_text(size = title_size, face = "bold"),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_text_size))
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
getCellSelectionPlotTitle <- function(cell_features_result, strTitlePrefix = "",
    transcriptFeature = "num_transcripts") {
    selected <- cell_features_result[which(cell_features_result$is_cell ==
        TRUE & cell_features_result[[transcriptFeature]] > 0), ]

    numSTAMPs <- dim(selected)[1]
    readsPerUMI <- NA

    minNumUmis <- min(selected[[transcriptFeature]])
    medianUMI <- round(median(selected[[transcriptFeature]]))
    minIntronic <- min(selected$pct_intronic)
    comma_formatter <- scales::label_comma()

    strTitle <- sprintf("%s, intronic>=%.2f\n%s Nuclei, %d+ UMIs, medUMIs %s",
        strTitlePrefix, minIntronic, comma_formatter(numSTAMPs), minNumUmis,
        comma_formatter(medianUMI))

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
plotExpressionVsIntronic <- function(cell_features, densityCenters = NULL,
    intronicThreshold = NULL, title = "All cell barcodes", point_size = 0.25,
    alpha = 0.25, useCellBenderFeatures = TRUE) {
    if (!useCellBenderFeatures) {
        p <- ggplot() + theme_void()  # empty plot
        return(p)
    }

    # TO MAKE R CMD CHECK HAPPY
    x <- y <- num_transcripts <- pct_intronic <- frac_contamination <- NULL

    p <- ggplot(cell_features, aes(x = log10(num_transcripts),
        y = pct_intronic, color = frac_contamination)) +
        ggrastr::rasterize(geom_point(size = point_size, alpha=alpha),dpi=900) +
        scale_color_gradientn(
            colors = c("light blue", "blue", "red"), values = c(0,0.5,1)) +
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
            geom_point(data = data.frame(x = log10(densityCenters$empty[1]),
                y = densityCenters$empty[2]),
                aes(x = x, y = y), color = "red", size = 3) +
            geom_point(data = data.frame(x = log10(densityCenters$cell[1]),
                y = densityCenters$cell[2]),
                aes(x = x, y = y), color = "red", size = 3)
    }

    if (!is.null(intronicThreshold)) {
        p <- p +
            geom_hline(yintercept = intronicThreshold, linetype = "dotted",
                color = "red", linewidth = 1)
    }
    return(p)
}


#' Plot the expression vs intronic content of the cells
#' @param cell_features A data frame containing the cell features.
#' @param bounds_empty A list containing the bounds for the empty cells.
#' @param bounds_non_empty A list containing the bounds for the non-empty cells.
#' @param strTitleOverride The title of the plot.  Overrides a more generic
#'   title.
#' @param cex.axis The size of the axis text.
#' @param cex.lab The size of the axis labels.
#' @param cex.main The size of the main title.
#' @importFrom graphics par smoothScatter title rect
#' @noRd
plotCellTypeIntervals <- function(cell_features, bounds_empty, bounds_non_empty,
    strTitleOverride = NULL, cex.axis = 0.6, cex.lab = 0.7, cex.main = 0.65) {
    par(mar = c(2, 2, 2, 1), mgp = c(0.8, 0.25, 0), tck = -0.02)

    # how many cell barcodes in the training set?
    num_training_cells <- length(which(cell_features$training_label_is_cell))

    smoothScatter(log10(cell_features$num_transcripts),
        cell_features$pct_intronic, xlim = log10_UMI_AXIS_RANGE_NEW,
        ylim = c(0, 1), xlab = "log10( UMIs )", ylab = "intronic",
        main = "", axes = TRUE, cex.axis = cex.axis, cex.lab = cex.lab)

    strTitle <- paste("SVM Training -- initial nuclei exemplars selected [",
        num_training_cells, "]", sep = "")
    # optionally override title.
    if (!is.null(strTitleOverride)) {
        strTitle <- strTitleOverride
    }

    title(main = strTitle, line = 1, col.main = "black", cex.main = cex.main)

    graphics::rect(bounds_empty$umi_lower_bound,
        bounds_empty$intronic_lower_bound,
        bounds_empty$umi_upper_bound,
        bounds_empty$intronic_upper_bound,
        border = "red", lwd = 3, lty = 1)

    graphics::rect(bounds_non_empty$umi_lower_bound,
        bounds_non_empty$intronic_lower_bound,
        bounds_non_empty$umi_upper_bound,
        bounds_non_empty$intronic_upper_bound,
        border = "green", lwd = 3, lty = 1)

    # reset par to defaults.  trying to save the original PAR and
    # reset it was breaking.
    par(mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0), tck = NA)
}

plotSelectedCells <- function(cell_features_result, size = 0.25, alpha = 0.25) {
    strTitle <- "Selected Nuclei"

    df <- cell_features_result[cell_features_result$is_cell ==TRUE,]
    umi_min_threshold <- min(log10(df[["num_transcripts"]]))
    intronic_min_threshold <- min(df[["pct_intronic"]])

    # TO MAKE R CMD CHECK HAPPY
    num_transcripts <- pct_intronic <- is_cell <- NULL

    p <- ggplot(cell_features_result, aes(x = log10(num_transcripts),
        y = pct_intronic, color = is_cell)) +
        ggrastr::rasterize(geom_point(size = size, alpha=alpha),dpi=900) +
        labs(
            x = "log10(UMI)",
            y = "% Intronic",
            color = "Selected Cell"
        ) +
        ggtitle(strTitle) +
        scale_color_manual(values =
            c("TRUE" = "green", "FALSE" = "lightblue")) +
        coord_cartesian(xlim = log10_UMI_AXIS_RANGE_NEW) +
        theme_minimal() +
        guides(color = guide_legend(override.aes = list(size = 4, alpha=1)))

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
plotSelectedCellsSmoothScatter <- function(cell_features, strTitlePrefix = "",
    transcriptFeature = "num_transcripts", changePar = TRUE,
    useCellBenderFeatures = TRUE) {

    if (!useCellBenderFeatures) {
        p <- ggplot2::ggplot() + ggplot2::theme_void()  # empty plot
        return(p)
    }

    xlab <- "log10 ( UMIs )"
    strPrefix <- ""

    if (transcriptFeature == "num_retained_transcripts") {
        xlab <- "log10 ( UMIs post CBRB )"
    }

    cell_features_filtered <- cell_features[!is.na(cell_features$is_cell),
        ]
    # for cases where plotting remove background processed data.
    # number of transcripts should never be 0.
    cell_features_filtered <-
        cell_features_filtered[cell_features_filtered[[transcriptFeature]] >0, ]

    strTitle <- getCellSelectionPlotTitle(cell_features_filtered,
        strTitlePrefix = strTitlePrefix, transcriptFeature = transcriptFeature)
    df <- cell_features_filtered[cell_features_filtered$is_cell ==TRUE &
                                     cell_features_filtered[[transcriptFeature]] > 0, ]
    umi_min_threshold <- min(df[[transcriptFeature]])
    intronic_min_threshold <- min(df[["pct_intronic"]])

    if (changePar) {
        opar <- graphics::par(no.readonly = TRUE)
        graphics::par(mar = c(2, 2, 2, 1), mgp = c(0.8, 0.25, 0), tck = -0.02)
    }

    graphics::smoothScatter(log10(cell_features[[transcriptFeature]]),
        cell_features$pct_intronic, xlim = log10_UMI_AXIS_RANGE_NEW, ylim = c(0,
            1), xlab = xlab, ylab = "intronic", axes = TRUE, cex.axis = 0.6,
        cex.lab = 0.7)

    graphics::title(main = strTitle, line = 0.25,
        col.main = "black", cex.main = 0.65)
    graphics::abline(v = log10(umi_min_threshold), col = "red", lty = 2)
    graphics::abline(h = intronic_min_threshold, col = "red", lty = 2)

    if (changePar) {
        graphics::par(opar)
    }
}

plotScaledTrainingDataFeatures <- function(trainingData) {
    # This is to avoid using dplyr and adding more dependencies.
    # Find the index of the label column
    label_col_index <- which(names(trainingData) == "training_label_is_cell")

    # Exclude the label column when generating feature names and
    # values
    feature_columns <- trainingData[, -label_col_index]
    label_column <- trainingData[, label_col_index]

    # Create the long-format data frame
    longData <- data.frame(feature = rep(names(feature_columns),
        each = nrow(trainingData)), value = unlist(feature_columns),
        is_nuclei = rep(label_column, times = ncol(feature_columns)))

    # Plot TO MAKE R CMD CHECK HAPPY
    feature <- value <- is_nuclei <- NULL

    p<- ggplot(longData, aes(x = feature, y = value, fill = is_nuclei)) +
        geom_violin(trim = FALSE, alpha = 0.7) +
        scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
        labs(title = "Violin Plot of Training Data Features",
            x = "Features",
            y = "Values",
            fill = "Is Nuclei") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "top") +
        coord_flip()

    return(p)
}


plotCellProbabilities <- function(cell_features,
    strTitle = "Nuclei Probability") {

    df <- cell_features[which(cell_features$is_cell_prob >= 0.5), ]

    breaks <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1)
    colors <- c("red", "purple", "orange", "blue", "green")

    # Assign bins to is_cell_prob
    df$is_cell_prob_bin <-
        cut(df$is_cell_prob, breaks = breaks, include.lowest = TRUE)

    # Plot

    # TO MAKE R CMD CHECK HAPPY
    num_transcripts <- pct_intronic <- is_cell_prob_bin <- NULL

    p <- ggplot(df, aes(x = log10(num_transcripts), y = pct_intronic,
        color = is_cell_prob_bin)) +
        ggrastr::rasterize(geom_point(size = 0.5, alpha=0.25), dpi=900) +
        labs(
            x = "log10(UMI)",
            y = "% Intronic",
            color = "Cell Probability"
        ) +
        #coord_cartesian(xlim = log10_UMI_AXIS_RANGE_NEW) +
        ggtitle(strTitle) +
        theme_minimal() +
        scale_color_manual(values = colors) +
        guides(color = guide_legend(override.aes = list(size = 5, alpha=1)))
    return(p)
}



####################### GENE MODULE PLOTS

scatterPlotModuleScore <- function(cell_features,
    moduleName = "nuclei_gene_module_score", strTitle = "",
    point_size = 0.25, alpha = 0.25) {

    # TO MAKE R CMD CHECK HAPPY
    num_transcripts <- pct_intronic <- NULL

    if (!moduleName %in% colnames(cell_features)) {
        stop("Module name ", moduleName, " not found in cell features")
    }

    p <- ggplot(cell_features, aes(x = log10(num_transcripts),
        y = pct_intronic, color = get(moduleName))) +
        ggrastr::rasterize(geom_point(
            size = point_size, alpha = alpha), dpi = 900) +
        scale_color_gradient2(
            low = "red", mid = "lightgrey", high = "blue", midpoint = 0) +
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

plotModuleVsFracContam <- function(cell_features,
    moduleName = "nuclei_gene_module_score", point_size = 0.5, alpha = 0.5) {
    # TO MAKE R CMD CHECK HAPPY
    frac_contamination <- pct_intronic <- NULL

    if (!moduleName %in% colnames(cell_features)) {
        stop("Module name ", moduleName, " not found in cell features")
    }

    p <- ggplot(cell_features, aes(x = frac_contamination, y=get(moduleName),
        color = pct_intronic)) +
        ggrastr::rasterize(geom_point(size = point_size, alpha = alpha),
        dpi=900) +
        labs(
            x = "Fraction Contamination",
            y = moduleName,
            color = "% Intronic"
        ) +
        theme_minimal()
    return(p)
}

plotCellProbabilityConditionalCbrb <- function(cell_features_result,
    frac_contamination_threshold = 1,
    useCellBenderFeatures = TRUE) {

    # if not using celbender, short circuit with an empty plot
    if (!useCellBenderFeatures) {
        return(ggplot() + theme_minimal())
    }

    df <- cell_features_result[
        cell_features_result$frac_contamination >= frac_contamination_threshold,
    ]

    # Filter out NAs
    df <- df[!is.na(df$is_cell_prob), ]

    # Calculate the number of cell barcodes
    numCellBarcodes <- dim(df[df$is_cell == TRUE, ])[1]

    r <- graphics::hist(df$is_cell_prob, breaks = seq(0, 1, by = 0.01),
        plot = FALSE)

    dd <- data.frame(countLog10 = log10(r$counts + 1), midpoint = r$mids)

    # TO MAKE R CMD CHECK HAPPY
    midpoint <- countLog10 <- NULL

    title<-paste("SVM Probability for CBRB empty cell barcodes\n[",
        numCellBarcodes, "] barcode selected as nuclei")

    p <- ggplot(dd, aes(x = midpoint, y = countLog10)) +
        geom_bar(stat = "identity", width = 0.01) +
        geom_vline(xintercept = 0.5, color = "red",
            linetype = "dashed", linewidth = 1) +
        labs(x = "Cell Probability", y = "Log10(Count + 1)") +
        ggtitle(title) +
        theme_minimal() +
        annotation_logticks(sides = "l")
    return(p)
}
