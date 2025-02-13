# Some global parameters
MIN_UMIs_PER_STAMP <- 20  # what's worth even including in analysis
log10_UMI_AXIS_RANGE_NEW <- c(log10(MIN_UMIs_PER_STAMP), 6)   # for plotting
random.seed <- 1

#' Call nuclei using an SVM trained on cell summary features.
#'
#' A convenience method for calling cells using an SVM trained on cell summary features.
#'
#' @param datasetName The name of the dataset.  Used for plotting.  Useful to use a unique experimental identifier (UEI) if available.
#' @param cellFeaturesFile The cell features file.  This file must contain a column cell_barcode that is used to match the cell barcodes in the DGE matrix
#' and a num_reads column that defines the number of uniquely mapped reads for the cell.  Additional expected columns
#' used to train the SVM are defined by the features vector.
#' @param dgeMatrixFile A file containing a dense DGE matrix in DropSeq format, or the directory where files in the MTX format exist.
#' If a directory is supplied, 3 files are expect in the directory: matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz that encode
#' the counts in sparse format, the genes, and the cell barcodes respectively.
#' @param optimusH5File The h5ad file produced by the Optimus pipeline.  This input contains both the DGE matrix and the cell features required for cell calling.  Cell level metrics
#' are appropriately transformed to match SVM expectations for training.
#' @param cellProbabilityThreshold The probability threshold for selecting cells.  This value scales from 0-1, with 1 being extremely confident/stringent.  Set to null to use the SVM defaults (>0.5 = nuclei).
#' This modifies plotting outputs and the output features file classification, but does not filter any cells from the output.
#' @param max_umis_empty All cell barcodes with fewer than this number of UMIs are not considered for any aspects of nuclei selection.  It is assumed
#' these cell barcodes have fewer UMIs than the ambient RNA.
#' @param features A list of features to use for cell selection.  The features are used to train the SVM.  These features are columns
#' of the cell features data frame, with the exception of empty_gene_module_score, which is generated from the gene expression data.
#' @param useCBRBFeatures When true, the cell bender feature frac_contamination is used for cell selection.
#' When false, these features are not used.  This modifies the features argument.
#' @param forceTwoClusterSolution When true, the initialization of the SVM will attempt to find a solution with two clusters.
#' In cases where an experiment is overloaded with nuclei, this may correct the initial set of nuclei and empty droplets selected.
#' This argument is specific to the density-style selection of exemplars that is only applicable when useCBRBFeatures is false.
#' @param outPDF The PDF file to write the plots to.
#' @param outFeaturesFile The cell features dataframe, further annotated by the SVM to include the cell probability and label,
#' along with which cell barcodes were used for training.
#' @param outCellBenderInitialParameters If useCBRBFeatures is false and this parameter is not null, the SVM
#' is trained without cellbenber features, and the parameters expected_cells and total_droplets_included are
#' estimated and written to this file.
#' @param random.seed The random seed to use for reproducibility.  Pass NA to disable setting the random seed.
#' @import logger
#' @importFrom utils write.table
#' @export
runIntronicSVM<-function (datasetName, cellFeaturesFile=NULL, dgeMatrixFile=NULL,
                          optimusH5File=NULL,
                          cellProbabilityThreshold=NULL, max_umis_empty=50,
                          features=NULL,
                          useCBRBFeatures=TRUE, forceTwoClusterSolution=FALSE,
                          outPDF=NULL, outFeaturesFile=NULL, outCellBenderInitialParameters=NULL, random.seed=1) {

    if (!is.na(random.seed)) {
        set.seed(random.seed)
    }
    #parse and validate inputs.
    r=parseInputs(cellFeaturesFile=cellFeaturesFile, dgeMatrixFile=dgeMatrixFile, optimusH5File=optimusH5File)
    cell_features=r$cell_features
    dgeMatrix=r$dgeMatrix
    svmNucleusCaller = SvmNucleusCaller(
      cellFeatures = cell_features,
      dgeMatrix = dgeMatrix,
      cellProbabilityThreshold = cellProbabilityThreshold,
      maxUmisEmpty = max_umis_empty,
      featureColumns = features,
      forceTwoClusterSolution = forceTwoClusterSolution,
      useCBRBFeatures = useCBRBFeatures,
        datasetName = datasetName
    )

    cell_features_result=svmNucleusCaller$cell_features

    # Open the PDF device
    if (!is.null(outPDF)) {
        #pdf(outPDF, width = pdf_width, height = pdf_height)
        grDevices::pdf(outPDF)
    }

    if (!is.null(outPDF)) {
        grDevices::dev.off()
    }

    if (!is.null(outFeaturesFile)) {
        write.table(cell_features_result, outFeaturesFile, row.names=F, col.names=T, quote=F, sep="\t")
    }

    #get the initial cellbender remove-background parameters and write to a file.
    if (!is.null(outCellBenderInitialParameters)) {
        df=cell_features_result
        threshold_total_droplets=round(mean(df[df$training_label_is_cell==F,]$num_transcripts, na.rm=T))
        result=data.frame(total_droplets_included=length(which(df$num_transcripts>threshold_total_droplets)), expected_cells=length(which(df$is_cell==T)))
        write.table(result, outCellBenderInitialParameters, row.names=F, col.names=T, quote=F, sep="\t")
    }
}



#' Select cells using an SVM trained on cell summary features
#'
#' @param dataset_name The name of the dataset.  Used for plotting.  Useful to use a unique experimental identifier (UEI) if available.
#' @param cell_features A dataframe of cell features fproduced by buildCellFeaturesSimple()
#' @param dgeMatrix A dense or sparse DGE matrix.
#' @inheritParams runIntronicSVM
#' @return A list containing the dataset name, the cell features with the SVM results, various plots, and the DGE matrix if it was parsed.
#' @import e1071 hdrcde ggplot2 grid gridExtra cowplot ggrastr logger gridGraphics
#' @noRd
callByIntronicSVM<-function (dataset_name, cell_features, dgeMatrix,
                             cellProbabilityThreshold=NULL, max_umis_empty=50,
                             features=NULL,
                             forceTwoClusterSolution=FALSE) {

    #validate the features are in the cell_features data frame.
    #do not look at empty_gene_module_score, as it is constructed.
    for (f in setdiff(features, "empty_gene_module_score")) {
        if (!(f %in% colnames(cell_features))) {
            strLog=paste("Requested Features [", paste (features, collapse=", "), "]. Feature", f, "not found in cell_features data frame.")
            stop(strLog)
        }
    }

    #CBRB specific parameter, hard-coded, used by multiple functions if CBRB mode is enabled.
    maxContaminationThreshold=0.1

    log_info("Beginning SVM Intronic Cell Selection")
    log_info(paste("Features used for cell selection", paste(features, collapse=", ")))

    #swap the cell barcodes into the row names.
    rownames(cell_features) <- cell_features$cell_barcode
    cell_features$cell_barcode <- NULL

    #check if we are using the CBRB features
    useCellBenderFeatures="frac_contamination" %in% features

    #Gather the training data bounds
    allBounds=findTrainingDataBounds(cell_features, max_umis_empty=max_umis_empty, useCellBenderFeatures=useCellBenderFeatures, forceTwoClusterSolution=forceTwoClusterSolution)
    bounds_empty=allBounds$bounds_empty
    bounds_non_empty=allBounds$bounds_non_empty

    #label cell vs not for training.  Test data is NA.  Bounds are always "extended" so get rid of the extra plot annotations.
    cell_features_labeled=labelTrainingData(cell_features, bounds_empty, bounds_non_empty,
                                            maxContaminationThreshold=maxContaminationThreshold, useCellBenderFeatures=useCellBenderFeatures)
    logTrainingDataSelection(cell_features_labeled)

    #Add the gene module score to the cell features.
    #generate the gene module plots.
    r=addGeneModules (cell_features_labeled=cell_features_labeled, dgeMatrix=dgeMatrix, numGenes=100, useCellBenderFeatures=useCellBenderFeatures, verbose=F)
    cell_features_labeled=r$cell_features
    geneModulePlots=r$gene_module_plots[1:4]

    #in the rare case where there are no differentially expressed genes to generate a module score
    #drop the feature from downstream analysis.
    if (!"empty_gene_module_score" %in% colnames (cell_features_labeled)) {
        log_warn("empty_gene_module_score was not computed. Dropping feature from SVM analysis.")
        features=setdiff(features, "empty_gene_module_score")
    }

    #modal - if not using the cellbender features, this plot is empty.
    p1 = plotExpressionVsIntronic(cell_features_labeled, title = "All cell barcodes", useCellBenderFeatures=useCellBenderFeatures)

    p2 = cowplot::ggdraw(function() plotCellTypeIntervals(cell_features_labeled, bounds_empty, bounds_non_empty))

    # Run the SVM - this scales the feautures, trains the model, and predicts all cell barcode labels/probabilties.
    svm_result <- runSVM (cell_features_labeled, features, bounds_empty=bounds_empty, cellProbabilityThreshold=cellProbabilityThreshold)

    cell_features_result=svm_result$cell_features
    trainingData=svm_result$trainingData

    #plotting
    p3=plotSelectedCells(cell_features_result)

    strTitlePrefix = "SVM nuclei method"
    #add the ambient peak to the title
    ambientPeak=round(median (cell_features_result[which(cell_features_result$is_cell==F),]$num_transcripts))
    strAmbientPeak = sprintf(", ambient peak(UMIs) %d", ambientPeak)

    p4 = plotCellProbabilities (cell_features_result, strTitle="Cell Probability")
    p5 <- ggdraw(function() plotSelectedCellsSmoothScatter(cell_features_result, transcriptFeature="num_transcripts",
                                                           strTitlePrefix = paste0(strTitlePrefix, strAmbientPeak)))

    p6=ggplot() + theme_void() #empty plot
    if (useCellBenderFeatures) {
        p6 <- ggdraw(function() plotSelectedCellsSmoothScatter(cell_features_result, transcriptFeature="num_retained_transcripts",
                                                               strTitlePrefix = strTitlePrefix, useCellBenderFeatures=useCellBenderFeatures))
    }

    #this new plot is only tracked for the DGE enabled SVM
    #unforunately it needs the full SVM to be run before it can be calculated
    if (!is.null(geneModulePlots)) {
        geneModulePlots[["cell_probability_histogram"]]=plotCellProbabilityConditionalCbrb(cell_features_result, frac_contamination_threshold=1, useCellBenderFeatures=useCellBenderFeatures)
    }

    feature_plot=plotScaledTrainingDataFeatures(trainingData[,c(features, "training_label_is_cell")])
    log_info("Nuclei selection finished")
    plot_list=list(cellbender=p1, initialization=p2, selected_nuceli=p3,nuceli_probability=p4,selected_nuclei_density=p5,selected_nuclei_density_rb=p6)
    return (list(dataset_name=dataset_name, cell_features=cell_features_result, features=features,
                 plots=plot_list, featurePlot=feature_plot,
                 geneModulePlots=geneModulePlots, bounds_empty=bounds_empty, bounds_non_empty=bounds_non_empty))
}

##################
# INPUT VALIDATION
##################

parseInputs<-function (cellFeaturesFile=NULL, dgeMatrixFile=NULL, optimusH5File=NULL) {
    if (!is.null(cellFeaturesFile) & !is.null (dgeMatrixFile) & !is.null(optimusH5File)) {
        stop("Supply either the cellFeaturesFile and dgeMatrixFile, or optimusH5File.  Supplying all 3 is ambiguous.")
    }

    if (is.null(cellFeaturesFile) & is.null (dgeMatrixFile) & is.null(optimusH5File)) {
        stop("Supply either the cellFeaturesFile and dgeMatrixFile, or optimusH5File.")
    }

    if (!is.null(optimusH5File)) {
        r=parseOptimusH5ad(optimusH5File)
        result=list(cell_features=r$cell_features, dgeMatrix=r$dge)
        return (result)
    }


    cell_features <- readCellFeatures(cellFeaturesFile)
    dgeMatrix=readDgeFile (dgeMatrixFile, cell_features)
    result=list(cell_features=cell_features, dgeMatrix=dgeMatrix)
    return (result)
}

###############################################
# FIND EXEMPLAR CLASS BOUNDS
###############################################

#' Find exemplars for nuclei and empty droplets
#' @param cell_features The cell features data frame
#' @param max_umis_empty Exclude cell barcodes from analysis with fewer than these many UMIs.  This threshold should
#' exclude noisy barcodes, but not exclude the empty droplet cloud.
#' @param useCellBenderFeatures When true, the CellBender remove-background feature frac_contamination is used for cell selection.
#' @param forceTwoClusterSolution When true, the function will attempt to find a solution with two clusters.  This may be useful when
#' the data is overloaded and breaks the normal assumptions, but may find suboptimal solutions for other data sets.
#' @param verbose Print verbose output to log
#' @return A list containing the bounds for the empty droplet and nuclei exemplars
#' @import logger
#' @noRd
findTrainingDataBounds<-function (cell_features, max_umis_empty=50, useCellBenderFeatures=T, forceTwoClusterSolution=F, verbose=F) {
    if (useCellBenderFeatures) {
        return (findTraingDataBoundsCBRB(cell_features, max_umis_empty=max_umis_empty,))
    }

    #This is more complicated because it supports multiple models.
    max_umis_empty_off=20

    if (forceTwoClusterSolution) {
        #with the default UMI filter, looking for the separation between the two highest peaks
        log_info(paste("Exemplar selection [pit between peaks] with default UMI filter [", max_umis_empty, "]" ))
        defaultPitBetween=findTrainingDataBoundsDefaultIterative(cell_features, max_umis_empty=max_umis_empty, method="PitBetweenHighestPeaks", verbose=verbose)

        # #with the default UMI filter, looking for the separation between the two highest peaks
        log_info(paste("Exemplar selection [pit between peaks] with UMI filter [", max_umis_empty_off, "]" ))
        noUmiFiltertPitBetween=findTrainingDataBoundsDefaultIterative(cell_features, max_umis_empty=max_umis_empty_off, method="PitBetweenHighestPeaks", verbose=verbose)
        results=list("Separate Two Highest Peaks"=defaultPitBetween, "Separate Two Highest Peaks No Filter"=noUmiFiltertPitBetween)
        #In the forced two cluster solution, do not enforce some number of empty droplets.
        result=selectBestTrainingBoundsModel(results)
    } else {
        #the default UMI filter as set by the user.
        log_info(paste("Exemplar selection [pit after highest peak] with default UMI filter [", max_umis_empty, "]"))
        default=findTrainingDataBoundsDefaultIterative(cell_features, max_umis_empty=max_umis_empty, verbose=verbose)

        #without a UMI filter, for low noise experiments
        log_info(paste("Exemplar selection [pit after highest peak] with UMI filter [", max_umis_empty_off, "]"))
        noUmiFilter=findTrainingDataBoundsDefaultIterative(cell_features, max_umis_empty=max_umis_empty_off, verbose=verbose)
        results=list("Default Selection"=default, "No UMI Filter Selecton"=noUmiFilter)
        #If not forcing a two cluster model, enforce some minimum fraction of empty droplets!
        result=selectBestTrainingBoundsModel(results, minEmptyDropletFraction=0.2)
    }


    return (result)
}

#' Select the best model from a list
#'
#' This function selects the best model from a list of results.  If the minimum empty droplet fraction is set, it will
#' only consider models that have a minimum fraction of empty droplets.
#' If no models meet the minimum empty droplet fraction, it will select the best model from the original list.
#'
#' @param results A list of results from findTrainingDataBoundsDefaultIterative
#' @param minEmptyDropletFraction The minimum fraction of empty droplets that must be present in the training data.
#' @return The best model from the list
#' @import logger
#' @noRd
selectBestTrainingBoundsModel <- function(results, minEmptyDropletFraction=NULL) {
    # Extract the best_silhouette score from each result
    silhouettes <- sapply(results, function(res) res$best_silhouette)
    numEmpty <- sapply(results, function(res) res$numEmpty)
    numNonEmpty<- sapply(results, function(res) res$numNonEmpty)

    #the row names of the data frame are the model names.
    df=data.frame(silhouettes=silhouettes, numEmpty=numEmpty, numNonEmpty=numNonEmpty)

    #remove NA silhouette results.
    # this explicitly ignores NAs, in cases where the model failed to find good bounds.
    df=df[!is.na(df$silhouettes),]

    #get the best silhouette score name
    best_name <- rownames(df[which.max(df$silhouettes),])

    #if maxEmptyDropletFraction is not null, drop any results that have a low number of empty droplets.
    if (!is.null(minEmptyDropletFraction)) {
        df=df[df$numEmpty/df$numNonEmpty>minEmptyDropletFraction,]
    }

    #if there are results left with reasonable numbers of empty/non empties, keep the best silhouette score.
    #if there are no results left, keep the best silhouette score from the original list.
    if (nrow(df)>0) {
        best_name <- rownames(df[which.max(df$silhouettes),])
    }

    # Log the selected result
    log_info(paste("Selected result [", best_name,
                   "] with best silhouette score of [", round(max(silhouettes),3), "]", sep=""))

    #label the results that have the best initialization
    result=results[[best_name]]
    result$method=best_name

    # Return the best name and result as a list
    return(result)

}

findTraingDataBoundsCBRB<-function (cell_features, max_umis_empty=50, maxContaminationThreshold=0.1) {
    #the interval for empty cells for training data
    #explicitly avoid searching in the very low UMI area, which are likely tons of PCR error barcodes.
    cell_features_empty <- cell_features[cell_features$frac_contamination==1 & cell_features$num_transcripts>max_umis_empty,]
    bounds_empty=getHighestDensityIntervalsEnforcedSmoothing(cell_features_empty, pctDensity = 75)

    #the intervals for nuclei examples
    b=selectNucleiExemplarBounds(cell_features, maxContaminationThreshold, max_umis_empty=max_umis_empty, initialDensity=95, bounds_empty=bounds_empty, extendCellSelectionBounds=T)
    bounds_non_empty=b$bounds_non_empty
    bounds_non_empty_extended=b$bounds_non_empty_extended

    result=list(bounds_empty=bounds_empty, bounds_non_empty=bounds_non_empty_extended)
    return (result)
}


#' Find the training data bounds using only the \% intronic and num_transcripts features
#'
#' This function uses a turning point method to find the largest peak and a peak after that
#' to define a linear separation of the data into the set that likely contains empty droplets
#' and the part of the data that contains non-empty droplets. It then selects areas of maximum
#' density within each of the two partitions.
#'
#' @param cell_features A dataframe containing the dataset with at least columns num_transcripts
#'   and pct_intronic.
#' @param max_umis_empty Exclude cell barcodes from analysis with fewer than this many UMIs.
#'   This is used to exclude very small barcodes that are noise and may confuse the initial
#'   turning point threshold.
#' @param umiThresholdOverride A numeric value to override the UMI threshold. If NULL, the function
#'   will find the UMI threshold using the turning point method. If not NULL, use an initial
#'   estimate provided by another means.
#' @param verbose Logical. if TRUE, prints verbose output to the log.
#'
#' @return A list with the training data bounds.
#' @noRd
findTrainingDataBoundsDefault<-function (cell_features, max_umis_empty=50, umiThresholdOverride=NULL, verbose=T) {

    #drop very small empty cell barcodes
    df=cell_features[cell_features$num_transcripts>=max_umis_empty,]

    #could I split cell barcodes into two groups using the the low density region in-between the peaks?
    #This will be used to subtract the ambient peak and low UMI cell barcodes
    #using the gridsearch to find non-NA results, which makes this more robust.
    umiThreshold=umiThresholdOverride
    if (is.null(umiThreshold))
        umiThreshold <- PitAfterHighestPeakWithGridSearch(log10(df$num_transcripts))

    if (verbose)
        log_info("Initial UMI Threshold [", round(umiThreshold,3), "]")

    df2=df[log10(df$num_transcripts)<umiThreshold,]

    #get the highest density interval, this is normally the ambient peak.
    bounds_empty=getEmptyCellsByDensity (cell_features=df2, yAxisFeature="pct_intronic", pctDensity=75, showPlot=F, verbose=verbose)

    #if unable to select bounds from this umiThreshold of the data, exit out.
    if (is.null(bounds_empty)) {
        bounds_empty$umi_upper_bound=umiThreshold
        bounds_non_empty=data.frame(umi_lower_bound=NA, umi_upper_bound=NA, intronic_lower_bound=NA, intronic_upper_bound=NA)
        result=list(bounds_empty=bounds_empty, bounds_non_empty=bounds_non_empty)
        return (result)
    }

    #consider all cells to the right of the division.
    df2=df[log10(df$num_transcripts)>umiThreshold,]

    #early exit if there's no data left.
    if (nrow(df2)==0) {
        log_warn("No data left after selecting empty barcodes.")
        bounds_empty$umi_upper_bound=umiThreshold
        bounds_non_empty=data.frame(umi_lower_bound=NA, umi_upper_bound=NA, intronic_lower_bound=NA, intronic_upper_bound=NA)
        result=list(bounds_empty=bounds_empty, bounds_non_empty=bounds_non_empty)
        return (result)
    }

    #select the maximum %intronic and UMIs density from the remaining data.
    #Use a tighter bound, and extract the %intronic.
    bounds_non_empty_intronic=getHighestDensityIntervalsEnforcedSmoothing (df2, yAxisFeature="pct_intronic", pctDensity=75, maxPeaksExpected=1, showPlot=F)
    df3=df2[df2$pct_intronic>=bounds_non_empty_intronic$intronic_lower_bound & df2$pct_intronic<=bounds_non_empty_intronic$intronic_upper_bound,]
    #then rerun on that intronic slide to get the num UMIs.
    bounds_non_empty_transcripts=getHighestDensityIntervalsEnforcedSmoothing (df3, yAxisFeature="pct_intronic", pctDensity=75, showPlot=F)

    #Extend the upper bound UMI bound to the 99th quantile
    topUMI=as.numeric (stats::quantile(log10(df3$num_transcripts), 0.995))
    bounds_non_empty_transcripts$umi_upper_bound=topUMI

    #merge the two bounds into a single data frame.
    bounds_non_empty=data.frame(umi_lower_bound=bounds_non_empty_transcripts$umi_lower_bound, umi_upper_bound=bounds_non_empty_transcripts$umi_upper_bound,
                                intronic_lower_bound=bounds_non_empty_intronic$intronic_lower_bound, intronic_upper_bound=bounds_non_empty_intronic$intronic_upper_bound)

    result=list(bounds_empty=bounds_empty, bounds_non_empty=bounds_non_empty)
    return (result)

}

#' Find the training data bounds using an iterative approach
#'
#' This replaces findTrainingDataBoundsDefault with an iterative approach to find the training data bounds.
#' This function performs an iterative process to select the initial empty/nuclei clusters
#' by adjusting the density smoothing factor, computing UMI thresholds, and optimizing
#' based on silhouette scores. The function applies a grid search over the density adjustments
#' to identify the optimal UMI threshold for separating nuclei and empty droplets.
#'
#' @param cell_features A dataframe containing an entry for each cell barcode with columns num_transcripts and pct_intronic.
#' @param max_umis_empty Exclude cell barcodes from analysis with fewer than this many UMIs. This threshold should
#'        exclude noisy barcodes, but not exclude the empty droplet cloud.
#' @param method A character string specifying the method to use for selecting the UMI threshold. Options are:
#'        PitAfterHighestPeak (default) and PitBetweenHighestPeaks. The former selects the UMI threshold
#'        after the highest peak in the density plot, while the latter selects the UMI threshold between the two
#'        highest peaks in the density plot.
#' @param smoothingMultiple A numeric multiplier that controls the rate of smoothing adjustment.
#'        Default is 1.1.
#' @param max_iterations Maximum number of iterations for adjusting the smoothing factor. Default is 10.
#' @param early_termination_threshold Numeric value representing the percentage drop in silhouette score
#'        for early termination. Default is 0.05.
#' @param verbose Print verbose output to log
#' @return A list containing
#' A dataframe with columns smoothingMultipliers, umiThresholds, and silhouetteScores for each iteration.
#' The bounds found for the maximum silhouette score and the maximum silhouette score.
#'
#' @import logger
#' @noRd
findTrainingDataBoundsDefaultIterative<- function (cell_features, max_umis_empty=50,
                method = c("PitAfterHighestPeak", "PitBetweenHighestPeaks"),
                smoothingMultiple = 1.1,
                max_iterations = 10,
                early_termination_threshold = 0.05, verbose=T) {

    # Ensure the method is one of the allowed options
    method <- match.arg(method)

    if (verbose)
        log_info("Starting iterative smoothing with silhouette analysis to select initial empty/nuclei clusters using method [", method, "]")

    # Filter data with num_transcripts >= 50 (moved outside the loop)
    if (verbose)
        log_info("Filtering data with num_transcripts >= ", max_umis_empty, " for training data bounds.")

    df_filtered <- cell_features[cell_features$num_transcripts >= max_umis_empty, ]

    # Log10 transform num_transcripts
    x <- log10(df_filtered$num_transcripts + 1)

    # Initialize lists to store results within the loop
    results <- list(
        smoothingMultiplierList = c(),
        umiThresholdList = c(),
        silhouetteList = c()
    )

    # Helper function to accumulate results
    accumulate_results <- function(results, smoothing, umi, silhouette, log_message = NULL, verbose=T) {
        results$smoothingMultiplierList <- c(results$smoothingMultiplierList, smoothing)
        results$umiThresholdList <- c(results$umiThresholdList, umi)
        results$silhouetteList <- c(results$silhouetteList, silhouette)
        if (!is.null(log_message) & verbose) log_info(log_message)
        return(results)  # Return the updated results list
    }

    # Initialize smoothing multiplier
    smoothingMultiplier <- 1

    # Track the best score so far
    best_silhouette <- -Inf
    best_bounds=NULL
    best_umi_threshold=NULL

    # Iterate through the smoothing and silhouette calculations
    for (iter in 1:max_iterations) {

        # Compute density adjustment range
        denRange <- seq(0.25, 3, 0.25) ^ smoothingMultiplier

        # Perform density estimation
        den <- density(x, adjust = mean(denRange), n = 200)

        # Find the UMI threshold using the selected method.
        if (method == "PitAfterHighestPeak")
            umiThreshold <- PitAfterHighestPeakWithGridSearch(x, adjust = denRange)
        else if (method == "PitBetweenHighestPeaks")
            umiThreshold <- PitBetweenHighestPeaksWithGridSearch(x, adjust = denRange)
        else
            stop("Invalid method. Must be one of 'PitAfterHighestPeak' or 'PitBetweenHighestPeaks'")

        #if the umiThreshold is NA, the gridsarch failed.  Maybe a higher bandwidth will work...
        #the findTrainingDataBoundsDefault function will fail with an NA umi threshold, so don't even try.
        if (is.na(umiThreshold)) {
            results <- accumulate_results(results, smoothingMultiplier, umiThreshold, NA,
                                          "PitAfterHighestPeakWithGridSearch returned NA. Trying next iteration.", verbose = verbose)
            next
        }

        # Get the bounds using custom findTrainingDataBoundsDefault function
        bounds <- findTrainingDataBoundsDefault(df_filtered, max_umis_empty = max_umis_empty, umiThresholdOverride = umiThreshold, verbose = F)

        if (is.null(bounds)) {
            strLog=paste("Bounds Empty UMI [", )
            results <- accumulate_results(results, smoothingMultiplier, umiThreshold, NA)
            next
        }

        #label the data with clusters for the silhouette score
        cell_features_labeled=labelTrainingData(df_filtered, bounds$bounds_empty, bounds$bounds_non_empty,
                                                maxContaminationThreshold=NULL, useCellBenderFeatures=F, verbose=F)

        # Calculate silhouette score
        silhouetteResult=calculate_silhouette(cell_features_labeled, showPlot = FALSE, verbose=F)
        score = silhouetteResult$mean_silhouette

        if (is.na(score)) {
            results <- accumulate_results(results, smoothingMultiplier, umiThreshold, NA,
                                          "Silhouette score is NA. Trying next iteration.", verbose = verbose)
            next
        }

        if (score > best_silhouette) {
            best_silhouette <- score
            best_bounds=bounds
            best_umi_threshold=umiThreshold
        }

        results <- accumulate_results(results, smoothingMultiplier, umiThreshold, score,
                                      paste("Average silhouette score:",
                                            round(score, 3),
                                            "at UMI threshold", round(umiThreshold, 2)), verbose = verbose)

        # Check for early termination: if current score is X% lower than the best score
        if (score < (1 - early_termination_threshold) * best_silhouette) {
            if (verbose)
                log_info(sprintf("Early termination: Silhouette score dropped by more than %.1f%%", early_termination_threshold * 100))
            break
        }

        # Update the best silhouette score - now done above.
        # best_silhouette <- max(best_silhouette, score)

        # Update the smoothing multiplier for the next iteration
        smoothingMultiplier <- smoothingMultiplier * smoothingMultiple
    }

    # Return a dataframe of results after all iterations or early termination
    resultDF <- data.frame(smoothingMultipliers = results$smoothingMultiplierList,
                           umiThresholds = results$umiThresholdList,
                           silhouetteScores = results$silhouetteList)

    if (is.null(best_umi_threshold)) {
        log_info("Unable to find good initialization. Returning Empty Bounds")
        best_bounds=list(bounds_empty=data.frame(umi_lower_bound=NA, umi_upper_bound=NA, intronic_lower_bound=NA, intronic_upper_bound=NA),
                         bounds_non_empty=data.frame(umi_lower_bound=NA, umi_upper_bound=NA, intronic_lower_bound=NA, intronic_upper_bound=NA))
        result=list(best_silhouette=NA, best_umi_threshold=NA, bounds_empty=best_bounds$bounds_empty, bounds_non_empty=best_bounds$bounds_non_empty,
                    resultDF=resultDF, numEmpty=NA, numNonEmpty=NA)
        return(result)
    }

    log_info("Best silhouette score: [", round(best_silhouette,3),"] at UMI threshold [", round(best_umi_threshold,2), "]")

    result=list(best_silhouette=best_silhouette, best_umi_threshold=best_umi_threshold,
                bounds_empty=best_bounds$bounds_empty, bounds_non_empty=best_bounds$bounds_non_empty, resultDF=resultDF)

    #add the number of selected nuclei and empties for the best bounds.
    cell_features_labeled=labelTrainingData(df_filtered, result$bounds_empty, result$bounds_non_empty,
                                            maxContaminationThreshold=NULL, useCellBenderFeatures=F, verbose=F)

    result$numEmpty=sum(!is.na(cell_features_labeled$training_label_is_cell) & !cell_features_labeled$training_label_is_cell)
    result$numNonEmpty=sum(!is.na(cell_features_labeled$training_label_is_cell) & cell_features_labeled$training_label_is_cell)

    return(result)

}

#' Calculate the silhouette score for a clustering of cells.
#'
#' @param cell_features_labeled must include a column named "training_label_is_cell" that contains the cell/empty labels (TRUE for cell, FALSE for empty). NA entries are not included in the training data set.
#' @param downsampleRate The fraction of the data to use for silhouette calculation. Default is 0.1.
#' @param showPlot If TRUE, a silhouette plot will be displayed.
#' @param verbose If TRUE, the mean silhouette score will be printed to the log
#' @param seed The random seed to use for downsampling the data.
#' @return The mean silhouette score for the clustering.
#' @importFrom cluster silhouette
#' @importFrom scales rescale
#' @noRd
calculate_silhouette<-function (cell_features_labeled, downsampleRate=0.1, showPlot = FALSE, verbose=F, seed=123) {
    #restrict to the training data
    d=cell_features_labeled[!is.na(cell_features_labeled$training_label_is_cell),]

    # downsample to 10% of the original data.
    set.seed(seed)
    d <- d[sample(nrow(d), nrow(d)*downsampleRate),]

    # Ensure there are two clusters (TRUE and FALSE values)
    d$cluster_numeric <- as.numeric(factor(d$training_label_is_cell))
    if (length(unique(d$cluster_numeric)) < 2) {
        if (verbose) log_info("Not enough clusters to calculate silhouette.")
        result=list(mean_silhouette=NA, data=d)
        return (result)
    }

    # Extract the features used for clustering
    data_to_cluster <- d[, c("num_transcripts", "pct_intronic")]

    #rescale the feaures to 0-1 scale
    data_to_cluster <- data.frame(
        num_transcripts = scales::rescale(log10(data_to_cluster$num_transcripts+1)),
        pct_intronic = scales::rescale(data_to_cluster$pct_intronic)
    )

    # Compute the silhouette score based on the training_label_is_cell column
    silhouette_scores <- cluster::silhouette(d$cluster_numeric, stats::dist(data_to_cluster))

    # Optionally show the silhouette plot
    if (showPlot) {
        plot(silhouette_scores, main = "Silhouette plot")
    }

    data_to_cluster$silhouette_score=silhouette_scores[,3]
    data_to_cluster$cluster_numeric=d$cluster_numeric

    # Calculate the mean silhouette score for the clustering
    mean_silhouette <- mean(data_to_cluster$silhouette_score)

    if (verbose) {
        log_info(sprintf("Mean silhouette score: %.3f", mean_silhouette))
    }

    result=list(mean_silhouette=mean_silhouette, data=data_to_cluster)
    return(result)
}

#' Density bounds selection with enforced unimodal distribution via smoothing
#'
#' If the data is bimodal at the default bandwidth, the bandwidth is increased until a unimodal distribution is found.
#'
#' @param cell_features A data frame containing the cell features.
#' @param yAxisFeature A string specifying the feature to use for the y-axis.
#' @param pctDensity The density percentile to use for the bounds.
#' @param maxPeaksExpected The maximum number of peaks expected in the data.  The data will be iteratively smoothed until
#' this many peaks detected.
#' @param showPlot A boolean indicating whether to show the plot - this may display multiple plots along the search space.
#' @import hdrcde
#' @return A data frame containing the bounds for the UMI and intronic features.
#' @noRd
getHighestDensityIntervalsEnforcedSmoothing<-function (cell_features, yAxisFeature="pct_intronic", pctDensity=50, maxPeaksExpected=1, showPlot=F) {
    probList=unique(c(25,50,75,90,95,99,pctDensity))
    probList=probList[probList<=pctDensity]

    umiBounds=getBoundsByDensity(x=log10(cell_features$num_transcripts), probList=probList, pctDensity=pctDensity, maxPeaksExpected=maxPeaksExpected, showPlot=showPlot)$intervals
    yBounds=getBoundsByDensity(x=cell_features[[yAxisFeature]], probList=probList, pctDensity=pctDensity, maxPeaksExpected=maxPeaksExpected, showPlot=showPlot)$intervals
    df=data.frame(umi_lower_bound=umiBounds[1], umi_upper_bound=umiBounds[2], intronic_lower_bound=yBounds[1], intronic_upper_bound=yBounds[2])
    return (df)
}


#Need the maximum density under each interval!
#if expecting more than on peak, may need to break ties selecting the peak with the highest density.
#it may be useful to override the bandwidth with the iterative bandwidth used for the entire experiment,
#instead of coming up with a new estimate here.
getBoundsByDensity<-function (x, probList, pctDensity, maxPeaksExpected=1, showPlot=F) {
    pct=paste0(pctDensity, "%")
    unimodal=FALSE

    #the same bandwidth used by hdrcde::hdr.
    h = hdrcde::hdrbw(hdrcde::BoxCox(x, 1), mean(probList))

    #multiply bandwidth by 2 each iteration
    mult=1
    while (!unimodal) {
        h=h*mult
        zHDR=hdrcde::hdr(x, prob = probList, h=h)$hdr
        if (showPlot)
            hdrcde::hdr.den(x, prob=probList,h=h)

        result=zHDR[pct,]
        result=result[!is.na(result)]
        if (length(result)<=2*maxPeaksExpected) {
            resultDF=list(intervals=result, bw=h)
            return (resultDF)
        } else {
            mult=mult*2
        }
    }
}

#if there are multiple peaks, only keep the largest peak.
getEmptyCellsByDensity<-function (cell_features, yAxisFeature="pct_intronic", pctDensity=75, showPlot=F, verbose=F) {

    probList=unique(c(25,50,75,90,95,99,pctDensity))
    probList=probList[probList<=pctDensity]

    #do this by density.  It's fully expected this distribution is bimodal (empty/nuclei.)
    x=log10(cell_features$num_transcripts+1)
    zHDR=getBoundsByDensity(x, probList=probList, pctDensity=pctDensity, maxPeaksExpected=2, showPlot=showPlot)

    #which interval contains the majority of the data?
    # Step 1: Split zHDR into pairs of intervals
    zHDR_intervals <- matrix(zHDR$intervals, ncol = 2, byrow = TRUE)

    #copied from the internals of hdrcde::hdr.den
    #we use the same bandwidth as used in getBoundsByDensity.
    d=density(x, bw=zHDR$bw, n=1001)
    maxDensityX=d$x[which.max(d$y)]

    # Step 2: Calculate the fraction of scores in X within each interval
    containsInterval<-function (interval, maxDensityX) {
        return (maxDensityX>=interval[1] & maxDensityX<=interval[2])
    }

    # Step 3: Identify the interval with the highest density
    max_idx <- which(apply(zHDR_intervals, 1, containsInterval, maxDensityX))
    #this step can fail when the initial cutpoint is poor.
    if (length(max_idx)==0) {
        if (verbose)
            log_warn("Unable to select the empty cell barcode highest density interval.")
        return (NULL)
    }
    umiBounds <- zHDR_intervals[max_idx, ]

    #select data from this interval and find the %intronic.
    df=cell_features[log10(cell_features$num_transcripts)>=umiBounds[1] & log10(cell_features$num_transcripts)<=umiBounds[2],]
    pctDensity=95
    probList=unique(c(25,50,75,90,95,99,pctDensity))
    probList=probList[probList<=pctDensity]
    yBounds=getBoundsByDensity(df[[yAxisFeature]], probList, pctDensity, maxPeaksExpected=1, showPlot=showPlot)$intervals

    result=data.frame(umi_lower_bound=umiBounds[1], umi_upper_bound=umiBounds[2], intronic_lower_bound=yBounds[1], intronic_upper_bound=yBounds[2])
    return (result)
}


selectNucleiExemplarBounds<-function (cell_features, maxContaminationThreshold=0.1, max_umis_empty=50, initialDensity=95, bounds_empty=NULL, extendCellSelectionBounds=TRUE) {
    #the interval for cell barcodes that are not empty for training data.
    #set a threshold of the lowest 25% of contamination for cells that have contamination < 1.
    cell_features_non_empty <- cell_features[cell_features$frac_contamination<1,]
    contaminationThreshold=stats::quantile(cell_features_non_empty$frac_contamination, probs=seq(0,1,0.01))["25%"]
    cell_features_non_empty <- cell_features_non_empty[cell_features_non_empty$frac_contamination<=contaminationThreshold,]
    bounds_non_empty=getHighestDensityIntervalsEnforcedSmoothing(cell_features_non_empty, pctDensity=95)

    #if requested, try to extend the initial selection to a larger region.
    bounds_non_empty_extended=NULL
    if (extendCellSelectionBounds & !is.null(bounds_empty)) {
        idx=which(log10(cell_features$num_transcripts) > (bounds_empty$umi_lower_bound +1 ) &
                      cell_features$pct_intronic>=bounds_non_empty$intronic_lower_bound &
                      cell_features$pct_intronic<=bounds_non_empty$intronic_upper_bound &
                      cell_features$frac_contamination <= maxContaminationThreshold)
        bounds_non_empty_extended=getHighestDensityIntervalsEnforcedSmoothing(cell_features[idx,], pctDensity=95, showPlot = F)
        bounds_non_empty_extended=merge_bounds(bounds_non_empty, bounds_non_empty_extended)
    }
    return (list(bounds_non_empty=bounds_non_empty, bounds_non_empty_extended=bounds_non_empty_extended))

}



###################################
# USE BOUNDS TO LABEL TRAINING DATA
###################################

#' Label the training data based on the given bounds
#'
#' @param cell_features A data frame containing the cell features.
#' @param bounds_empty A data frame containing the bounds for the empty cells.
#' @param bounds_non_empty A data frame containing the bounds for the non-empty cells.
#' @param maxContaminationThreshold The maximum contamination threshold for non-empty cells (only when using cellbender features.)
#' @param useCellBenderFeatures A boolean indicating whether to use CellBender features.
#' @param verbose A boolean indicating whether to print log messages.
#' @return A data frame containing the training data with a new column `training_label_is_cell`.
#' @noRd
labelTrainingData<-function (cell_features, bounds_empty, bounds_non_empty,
                             maxContaminationThreshold=0.1, useCellBenderFeatures=T, verbose=T) {

    if (useCellBenderFeatures) {
        result=(labelTrainingDataCBRB(cell_features, bounds_empty, bounds_non_empty, maxContaminationThreshold=maxContaminationThreshold, verbose=verbose))
    } else {
        result=labelTrainingDataDefault(cell_features, bounds_empty, bounds_non_empty, verbose=verbose)
    }
    return (result)
}

logTrainingDataSelection<-function (cell_features_labeled) {
    numEmpty=sum(!is.na(cell_features_labeled$training_label_is_cell) & !cell_features_labeled$training_label_is_cell)
    numNonEmpty=sum(!is.na(cell_features_labeled$training_label_is_cell) & cell_features_labeled$training_label_is_cell)
    log_info("Number of empty exemplars: [", numEmpty, "]")
    log_info("Number of nuclei exemplars: [", numNonEmpty, "]")
}

labelTrainingDataCBRB<-function (cell_features, bounds_empty, bounds_non_empty, bounds_non_empty_extended=NULL,
                             maxContaminationThreshold=0.1, verbose=T) {
    # Define a function to find indices based on the given bounds
    find_indices <- function(cell_features, bounds, maxContaminationThreshold=1, minContaminationThreshold=0) {
        which(
            log10(cell_features$num_transcripts) >= bounds$umi_lower_bound &
                log10(cell_features$num_transcripts) <= bounds$umi_upper_bound &
                cell_features$pct_intronic >= bounds$intronic_lower_bound &
                cell_features$pct_intronic <= bounds$intronic_upper_bound &
                cell_features$frac_contamination <= maxContaminationThreshold &
                cell_features$frac_contamination >= minContaminationThreshold
        )
    }

    # Get the indices for empty and non-empty classes
    #empty classes have no filtering on max contamination - they are just empty - all UMIs always subtracted from the ambient.
    idxEmpty <- find_indices(cell_features, bounds_empty, 1, 1)
    #non empty classes have a max contamination threshold
    if (!is.null(bounds_non_empty_extended)) {
        idxNonEmpty <- find_indices(cell_features, bounds_non_empty_extended, maxContaminationThreshold, 0)
    } else {
        idxNonEmpty <- find_indices(cell_features, bounds_non_empty, maxContaminationThreshold, 0)
    }
    all=sort(union(idxEmpty, idxNonEmpty))

    # Assign classes based on the indices
    training_data=cell_features
    training_data$training_label_is_cell <- NA
    training_data$training_label_is_cell[idxNonEmpty] <- TRUE
    training_data$training_label_is_cell[idxEmpty] <- FALSE

    if (verbose)
        log_info("Training empty/nuclei cell barcodes selected [using CBRB]")
    return (training_data)
}

labelTrainingDataDefault<-function (cell_features, bounds_empty, bounds_non_empty, verbose=T) {
    # Define a function to find indices based on the given bounds
    find_indices <- function(cell_features, bounds) {
        which(
            log10(cell_features$num_transcripts) >= bounds$umi_lower_bound &
                log10(cell_features$num_transcripts) <= bounds$umi_upper_bound &
                cell_features$pct_intronic >= bounds$intronic_lower_bound &
                cell_features$pct_intronic <= bounds$intronic_upper_bound
        )
    }

    # Get the indices for empty and non-empty classes
    #empty classes have no filtering on max contamination - they are just empty - all UMIs always subtracted from the ambient.
    idxEmpty <- find_indices(cell_features, bounds_empty)
    #non empty classes have a max contamination threshold
    idxNonEmpty <- find_indices(cell_features, bounds_non_empty)
    all=sort(union(idxEmpty, idxNonEmpty))

    # Assign classes based on the indices
    training_data=cell_features
    training_data$training_label_is_cell <- NA
    training_data$training_label_is_cell[idxNonEmpty] <- TRUE
    training_data$training_label_is_cell[idxEmpty] <- FALSE

    if (verbose)
        log_info("Training empty/nuclei cell barcodes selected [using density-only method]")
    return (training_data)
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

######################################
# RUN SVM
######################################

runSVM<-function (cell_features_labeled, features, bounds_empty, cellProbabilityThreshold=NULL) {

    #scale the requested cell features.
    cell_features_scaled=scaleFeatures(cell_features_labeled, features)

    #get the training data - restrict data to the labeled training data for the specific features we are using.
    trainingData=cell_features_scaled[!is.na(cell_features_scaled$training_label_is_cell),c(features, "training_label_is_cell")]
    trainingData$training_label_is_cell=as.factor(trainingData$training_label_is_cell)

    # tune hyper parameters
    #parameters <- tune.svm(training_label_is_cell ~ ., data = trainingData, type = 'C-classification', kernel = 'radial', gamma = 10^(-2:1), cost = 10^(-2:1), probability=TRUE)
    #svm_model=parameters$best.model

    #train the SVM
    svm_model <- e1071::svm(training_label_is_cell ~ ., data = trainingData, type = 'C-classification', kernel = 'radial', probability=TRUE)

    #predict all points
    predictions <- stats::predict(svm_model, cell_features_scaled[,features], probability = TRUE)

    # Extract the probabilities
    probabilities <- attr(predictions, "probabilities")

    # Calculate confidence scores
    # confidence_scores <- apply(probabilities, 1, function(x) max(x) / sum(x))

    # Combine predictions and confidence into a dataframe
    results_df <- data.frame(
        is_cell = predictions,
        is_cell_prob= probabilities[,"TRUE"]
        #confidence = confidence_scores
    )

    #optionally override cell probability
    if (!is.null(cellProbabilityThreshold)) {
        results_df$is_cell=ifelse(results_df$is_cell_prob>=cellProbabilityThreshold, TRUE, FALSE)
    }

    # Exclude barcodes lower than the lower bound of UMIs from classification.
    # There's some noisy barcodes with few reads that don't have enough data to be reliably classified
    # but we know they are not true nuclei.
    idx=which(log10(cell_features_labeled$num_transcripts)<bounds_empty$umi_upper_bound)
    results_df$is_cell_prob[idx]=NA
    #results_df$confidence[idx]=NA
    results_df$is_cell[idx]=FALSE

    # merge the results_df with the original cell_features
    cell_features_result <- cbind(cell_features_labeled, results_df)
    result=list(cell_features=cell_features_result, trainingData=trainingData, svm_model=svm_model)
    return (result)
}

scaleFeatures<-function (cell_features_labeled, features) {
    #make an explicity copy
    cell_features_scaled=data.frame(cell_features_labeled)
    #log 10 the larger features
    cell_features_scaled$num_transcripts=log10(cell_features_scaled$num_transcripts)

    #This feature isn't used by default but is scaled in case someone wants to experiment.
    if ("num_reads" %in% cell_features_scaled)
        cell_features_scaled$num_reads=log10(cell_features_scaled$num_reads)

    #scale the features prior to SVN
    cell_features_scaled[,features]=scale(cell_features_scaled[,features])

    return (cell_features_scaled)
}



#############################################
# GENE MODULE DISCOVERY AND FEATURE BUILDING
#############################################

#' Find Nuclei/Empty gene modules
#'
#' Given the training nuclei and empty cell barcodes, find the gene modules that
#' are more differentially expressed between the groups, and score all cells with those features.
#'
#' @param cell_features_labeled The cell_features_labelel dataframe is the cell features file produced by
#' ,  with an additional indicator column "training_label_is_cell" that is TRUE for
#' nuclei and FALSE for empty.  This value is set to NA for all other cell barcodes.
#' @param dgeMatrix The matrix of gene expression.
#' @param numGenes The number of genes to select for each module.
#' @param useCellBenderFeatures If true, additional plots are generated to compare the empty module score to cellbender fraction of UMIs removed.
#' @param verbose A boolean indicating whether to print log messages.
#' @return A list containing the cell features with the gene module scores and QC plots.
#' If there are no differentially expressed genes, the empty_gene_module_score will be set to NA and the plots will be NULL.
#' @import Seurat presto Matrix methods
#' @noRd
addGeneModules<-function (cell_features_labeled, dgeMatrix, numGenes=100, useCellBenderFeatures=TRUE, verbose=F) {
    log_info("Adding gene module score(s)")
    #The annoying importing of methods is to make R CMD CHECK happy.
    dgeMatrix=methods::as(dgeMatrix, "CsparseMatrix")

    #garbage collection to get rid of dense matrix.
    gc()

    #make sure the dgeMatrix is ordered in the same way as the metadata before merging.
    dgeMatrix=dgeMatrix[,rownames(cell_features_labeled)]

    #create a seurat object of the cells, with expression and single cell metrics.
    seurat_object <- CreateSeuratObject(counts = dgeMatrix)
    seurat_object=add_cell_metadata (seurat_object, cell_features_labeled)

    #subset to the training cell barcodes.
    seurat_object_trainng=seurat_object[,which(!is.na(seurat_object$training_label_is_cell))]

    #pseudobulk to detect differentially expressed genes, standard normalization
    seurat_object_pseudobulked=pseudobulkEmpties(seurat_object_trainng, showPlot = FALSE, verbose=verbose)
    seurat_object_pseudobulked <- NormalizeData(seurat_object_pseudobulked, normalization.method = "LogNormalize", scale.factor = 10000, verbose=verbose)

    #get the differentially expressed genes in each direction.
    deWilcoxPB<-de_wilcox(seurat_object_pseudobulked, fdrThreshold=0.05, min.pct = 0.25, logfc.threshold = 1)

    #select up to numGenes from each direction.
    upGenes=deWilcoxPB[deWilcoxPB$avg_log2FC>0,]
    downGenes=deWilcoxPB[deWilcoxPB$avg_log2FC<0,]
    downGenes=downGenes[order(downGenes$avg_log2FC, decreasing = F),]
    geneListUp=head (rownames(upGenes), numGenes)
    geneListDown=head (rownames(downGenes), numGenes)

    #if there are no genes in the list, exit early.
    if (length(geneListDown)==0) {
        log_warn("No differentially expressed genes found for empty cell barcodes.")
        cell_features=cell_features_labeled
        result=list(cell_features=cell_features, downGenes=downGenes, gene_module_plots=list(training_data=NULL, frac_contamination=NULL, empty_gene_module_score=NULL, empty_gene_module_score_vs_contam=NULL))
        return (result)
    }


    #go back to the single cell data for scoring modules.
    #seurat_object_pseudobulked=AddModuleScore(seurat_object_pseudobulked, features = list(geneListUp), name = "nuclei_gene_module_score")
    seurat_object_pseudobulked=AddModuleScore(seurat_object_pseudobulked, features = list(geneListDown), name = "empty_gene_module_score", verbose=verbose)

    #go back to the single cell data for scoring modules.
    #this will add the score to all data, not just the training data.
    seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000, verbose=verbose)

    #seurat_object=AddModuleScore(seurat_object, features = list(geneListUp), name = "nuclei_gene_module_score")
    seurat_object=AddModuleScore(seurat_object, features = list(geneListDown), name = "empty_gene_module_score", verbose=verbose)

    #create a dataframe that includes the module scores
    #clean up the module names to remove the training "1".
    cell_features=seurat_object[[]]
    colnames(cell_features) =gsub(x = colnames(cell_features) , pattern = "1", replacement = "")
    dropCols=c("orig.ident", "nCount_RNA", "nFeature_RNA", "training_identity")
    cell_features = cell_features[,!(colnames(cell_features) %in% dropCols)]

    #QC plots
    #p1=VlnPlot(seurat_object_pseudobulked, features = c("nuclei_gene_module_score1", "empty_gene_module_score1"), group.by  = "training_identity")
    p1=VlnPlot(seurat_object_pseudobulked, features = c("empty_gene_module_score1"), group.by  = "training_identity")
    p2=ggplot() + theme_void() #empty plot
    p6=ggplot() + theme_void() #empty plot

    p4=scatterPlotModuleScore(cell_features, moduleName = "empty_gene_module_score", strTitle="Empty Score")

    #Only generate these plots if CBRB is in use.
    if (useCellBenderFeatures) {
        p2 = scatterPlotModuleScore(cell_features, moduleName = "frac_contamination", strTitle = "CellBender fraction contamination")
        p6 = plotModuleVsFracContam(cell_features, moduleName = "empty_gene_module_score")
    }
    result=list(cell_features=cell_features, downGenes=downGenes, gene_module_plots=list(training_data=p1, frac_contamination=p2, empty_gene_module_score=p4, empty_gene_module_score_vs_contam=p6))
    return (result)

}

pseudobulkEmpties<-function (seurat_object, showPlot=FALSE, verbose=FALSE) {

    # Separate target (nuclei) and non-target cells (empties)
    target_cells <- subset(seurat_object, subset = training_identity == "nuclei")
    non_target_cells <- subset(seurat_object, subset = training_identity == "empty")

    # Extract expression matrices
    # this is backwards compatable with Seurat V4
    target_expr <- GetAssayData(target_cells, assay = "RNA", layer = "counts")
    non_target_expr <- GetAssayData(non_target_cells, assay = "RNA", layer = "counts")

    # Pseudobulk non-target cells to match UMI count of target cells
    pseudobulked_data <- pseudobulk_to_match_umi_fast_optimized(target_expr, non_target_expr)
    # Remove columns with 0 total counts
    non_zero_cols <- colSums(pseudobulked_data) > 0
    filtered_pseudobulked_data <- pseudobulked_data[, non_zero_cols,drop=F]

    # Assign names to pseudobulked cells
    pseudobulked_cell_names <- paste0("EMPTY_", seq_len(ncol(filtered_pseudobulked_data)))
    colnames(filtered_pseudobulked_data) <- pseudobulked_cell_names

    # Combine target and pseudobulked non-target data
    combined_counts <- cbind(target_expr, filtered_pseudobulked_data)

    # Create a new Seurat object with the combined data
    combined_seurat <- CreateSeuratObject(counts = combined_counts)

    # Calculate metadata for pseudobulked cells
    pseudobulk_nCount_RNA <- colSums(filtered_pseudobulked_data)
    pseudobulk_nFeature_RNA <- colSums(filtered_pseudobulked_data > 0)
    pseudobulk_orig_ident <- rep("SeuratProject", length(pseudobulked_cell_names))

    # Create an empty data frame with the same columns as target_cells@meta.data
    empty_metadata <- target_cells@meta.data[0, ]
    empty_metadata <- empty_metadata[rep(1, length(pseudobulked_cell_names)), ]

    # Assign row names and fill in the metadata
    rownames(empty_metadata) <- pseudobulked_cell_names
    empty_metadata$orig.ident <- pseudobulk_orig_ident
    empty_metadata$nCount_RNA <- pseudobulk_nCount_RNA
    empty_metadata$nFeature_RNA <- pseudobulk_nFeature_RNA
    empty_metadata$training_identity <- "empty"

    # Combine the metadata
    combined_metadata <- rbind(target_cells@meta.data, empty_metadata)

    # Add the metadata to the combined Seurat object
    combined_seurat <- AddMetaData(combined_seurat, metadata = combined_metadata)

    if (showPlot) {
        z=combined_seurat@meta.data

        #TO MAKE R CMD CHECK HAPPY
        nCount_RNA=training_identity=NULL;

        p=ggplot(z, aes(x = nCount_RNA, fill = training_identity)) +
            geom_density(alpha = 0.5) +
            labs(title = "Density Plot of nCount_RNA by Training Identity",
                 x = "nCount_RNA",
                 y = "Density") +
            theme_minimal()
        print (p)

    }
    log_info("Finished Constructing pseudobulked empties")
    return(combined_seurat)
}

#' Pseudobulk empty droplets to match UMI count of target cells
#'
#' Function to pseudobulk non-target cells to match UMI count of target cells without converting to dense format
#' Optimized pseudobulk function with shuffling.  This function is faster than the original pseudobulk function,
#' but produces a very similar list of top gene results.
#'
#' @param target_expr The target data set (nuclei)
#' @param non_target_expr The non target data set to pseudobulk (empty droplets)
#' @param seed Random seed for reproducibility
#' @param verbose Print progress messages
#'
#' @return A sparse matrix with pseudobulked data
#' @importFrom Matrix colSums rowSums
#' @noRd
pseudobulk_to_match_umi_fast_optimized <- function(target_expr, non_target_expr, seed = 42, verbose = FALSE) {
    log_info("Constructing pseudobulked empty droplets for gene module detection")
    set.seed(seed)

    # Get UMI counts for target and non-target cells
    target_umis <- Matrix::colSums(target_expr)
    shuffled_target_umis <- sample(target_umis)
    non_target_umi_counts <- Matrix::colSums(non_target_expr)
    non_target_indices <- sample(length(non_target_umi_counts))
    non_target_umi_counts_shuffled <- non_target_umi_counts[non_target_indices]

    # Initialize the pseudobulk matrix
    pseudobulk_matrix <- Matrix::Matrix(0, nrow = nrow(non_target_expr), ncol = length(shuffled_target_umis), sparse = TRUE)
    rownames(pseudobulk_matrix) <- rownames(non_target_expr)

    # Compute cumulative sums
    cumsum_target_umis <- cumsum(shuffled_target_umis)
    cumsum_non_target <- cumsum(non_target_umi_counts_shuffled)

    # Initialize start and end indices
    start_indices <- integer(length(shuffled_target_umis))
    end_indices <- integer(length(shuffled_target_umis))
    prev_end_index <- 0

    for (i in seq_along(shuffled_target_umis)) {
        target_cumsum <- cumsum_target_umis[i]
        end_index <- findInterval(target_cumsum, cumsum_non_target) + 1

        if (end_index > length(cumsum_non_target)) {
            if (verbose) message("Not enough non-target UMIs to match target UMIs.")
            break
        }

        start_index <- prev_end_index + 1
        start_indices[i] <- start_index
        end_indices[i] <- end_index
        prev_end_index <- end_index
    }

    # Remove zero entries (in case we ran out of non-target UMIs)
    valid_indices <- which(start_indices > 0 & end_indices >= start_indices)
    for (i in valid_indices) {
        selected_cells_indices <- non_target_indices[start_indices[i]:end_indices[i]]
        pseudobulk_matrix[, i] <- Matrix::rowSums(non_target_expr[, selected_cells_indices, drop = FALSE])
        if (verbose && i %% 100 == 0) {
            message("Pseudobulked cell barcode [", i, "] of [", length(shuffled_target_umis), "]")
        }
    }

    return(pseudobulk_matrix)
}


de_wilcox<-function (seurat_object, fdrThreshold=1, min.pct=0.25, logfc.threshold = 0.25) {
    de_results <- FindMarkers(
        seurat_object,
        ident.1 = "nuclei",
        ident.2 = "empty",
        test.use = "wilcox",
        min.pct = min.pct, # only return genes that are detected in at least 25% of nuclei or empty cells
        #min.diff.pct = 0.25,        # only return genes that are detected in at least 25% of nuclei or empty cells
        logfc.threshold = logfc.threshold, # only return genes that have at least 0.25 log fold change difference
        group.by = "training_identity",
        only.pos=FALSE
    )
    de_results=de_results[order(de_results$avg_log2FC, decreasing = T),]
    de_results=de_results[de_results$p_val_adj <= fdrThreshold,]
    return (de_results)

}

add_cell_metadata<-function (seurat_object, cell_features) {
    cell_features$training_identity <- ifelse(cell_features$training_label_is_cell, "nuclei", "empty")
    seurat_object=AddMetaData(seurat_object, metadata = cell_features)
    Idents(seurat_object) <-cell_features$training_label_is_cell
    return (seurat_object)
}


########################
# PLOTTING CODE
########################


#' Create a single page of cell selection plots
#'
#' @param plots A list of ggplot2 plots generated by callByIntronicSVM
#' @param geneModulePlots (Optional) A list of ggplot2 plots that capture the expression feautres of the gene modules
#' @param featurePlot (Optional) A ggplot2 plot that captures the feature values for the training data.
#' @param dataset_name The name of the dataset
#' @param outPDF The PDF file to write the plots to
#' @param useOpenPDF If TRUE, do not open/close the outPDF.  If false, opens a new PDF device, writes the plots, and closes the device.
#' @import cowplot ggplot2
#' @noRd
arrangeSVMCellSelectionPlots <- function(plots, geneModulePlots=NULL, featurePlot=NULL, dataset_name, outPDF, useOpenPDF = FALSE) {
    #plots 1,3,4 are ggplot2.

    plots[[1]]=plots[[1]]+custom_theme()
    plots[[3]]=plots[[3]]+custom_theme()
    plots[[4]]=plots[[4]]+custom_theme()

    # Arrange the plots into a grid
    plot_grid <- plot_grid(plotlist = plots, ncol = 2, nrow = 3)

    # Add the title
    title <- ggdraw() +
        draw_label(dataset_name, fontface = 'bold', size = 12, hjust = 0.5)

    # Combine the title and the plot grid
    final_plot <- plot_grid(title, plot_grid, ncol = 1, rel_heights = c(0.05, 1))

    # Set PDF dimensions
    #pdf_width <- 12
    #pdf_height <- 12

    # Open the PDF device
    if (!useOpenPDF & !is.null(outPDF)) {
        #pdf(outPDF, width = pdf_width, height = pdf_height)
        grDevices::pdf(outPDF)
    }

    print(final_plot)

    plots[[4]]=plots[[4]]+custom_theme(title_size = 12, axis_title_size = 10, axis_text_size = 10, legend_title_size = 10, legend_text_size = 10)

    print (plots[[4]])

    if (!is.null(featurePlot)) {
        print(featurePlot)
    }

    if (!is.null(geneModulePlots)) {
        print(arrangeSVMGeneModulePlots(geneModulePlots, dataset_name))
    }

    if (!useOpenPDF & !is.null(outPDF)) {
        grDevices::dev.off()
    }
}

#' Create set of cell selection plots to evaluate cell selectionwhen CBRB features are not used.
#'
#' @param plots A list of ggplot2 plots generated by callByIntronicSVM
#' @param geneModulePlots (Optional) A list of ggplot2 plots that capture the expression feautres of the gene modules
#' @param featurePlot (Optional) A ggplot2 plot that captures the feature values for the training data.
#' @param dataset_name The name of the dataset
#' @param outPDF The PDF file to write the plots to
#' @param useOpenPDF If TRUE, do not open/close the outPDF.  If false, opens a new PDF device, writes the plots, and closes the device.
#' @import cowplot ggplot2
#' @noRd
arrangeSVMCellSelectionPlotsNoCBRB<-function (plots, geneModulePlots=NULL, featurePlot=NULL, dataset_name, outPDF, useOpenPDF = FALSE) {
    # re-order plots to include some plots that would otherwise be empty due to CBRB features not being used.
    # 1. initial selection exemplars (plots[[2]])
    # 2. empty score distribution based on exemplars (geneModulePlots[["training_data]])
    # 3. empty score distribution (geneModulePlots[[2]])
    # 4. Selected cells (plots[[3]])
    # 5. Selected cell probability (plots[[4]])
    # 6. SmoothScatter final selection [[plots[[5]]]]

    pList=list(plots[[2]], geneModulePlots[["training_data"]]+custom_theme(),
               geneModulePlots[["empty_gene_module_score"]]+custom_theme(),
            plots[[3]]+custom_theme(), plots[[4]]+custom_theme(), plots[[5]])

    # Arrange the plots into a grid
    plot_grid <- plot_grid(plotlist = pList, ncol = 2, nrow = 3)

    # Add the title
    title <- ggdraw() +
        draw_label(dataset_name, fontface = 'bold', size = 12, hjust = 0.5)

    # Combine the title and the plot grid
    final_plot <- plot_grid(title, plot_grid, ncol = 1, rel_heights = c(0.05, 1))

    # Open the PDF device
    if (!useOpenPDF & !is.null(outPDF)) {
        #pdf(outPDF, width = pdf_width, height = pdf_height)
        grDevices::pdf(outPDF)
    }

    print(final_plot)

    if (!is.null(featurePlot)) {
        print(featurePlot)
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
    #all plots are ggplot2 plots.

    plots_custom_theme <- lapply(plots, function(plot) plot + custom_theme())
    #https://github.com/satijalab/seurat/issues/3584
    plots_custom_theme[["training_data"]]= plots[["training_data"]] & custom_theme()

    #final plots to include:
    plotList=c("training_data", "empty_gene_module_score", "empty_gene_module_score_vs_contam", "cell_probability_histogram")
    plots_custom_theme=plots_custom_theme[plotList]

    # Arrange the plots into a grid
    nrow=ceiling(length(plots_custom_theme)/2)
    plot_grid <- plot_grid(plotlist = plots_custom_theme, ncol = 2, nrow = nrow)

    # Add the title
    title <- ggdraw() +
        draw_label(dataset_name, fontface = 'bold', size = 12, hjust = 0.5)

    # Combine the title and the plot grid
    final_plot <- plot_grid(title, plot_grid, ncol = 1, rel_heights = c(0.05, 1))
    return (final_plot)

}

# A custom theme for ggplot2 plots to make the text size appropriate when all plots are put together on a single page.
custom_theme <- function(title_size = 8, axis_title_size = 6, axis_text_size = 6, legend_title_size = 6, legend_text_size = 6) {
    theme(
        plot.title = element_text(size = title_size, face = "bold"),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.title = element_text(size=legend_title_size),
        legend.text = element_text(size=legend_text_size)
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
getCellSelectionPlotTitle<-function (cell_features_result, strTitlePrefix="", transcriptFeature="num_transcripts") {
    selected=cell_features_result[which(cell_features_result$is_cell==T),]

    numSTAMPs=dim(selected)[1]
    readsPerUMI=NA

    minNumUmis=min(selected[[transcriptFeature]])
    medianUMI=round (median(selected[[transcriptFeature]]))
    minIntronic=min(selected$pct_intronic)
    comma_formatter <- scales::label_comma()

    strTitle=sprintf("%s, intronic>=%.2f\n%s Nuclei, %d+ UMIs, medUMIs %s",
                     strTitlePrefix, minIntronic, comma_formatter(numSTAMPs), minNumUmis, comma_formatter(medianUMI))

    #if the number of reads is available, add the average reads/UMI to the title.
    if ("num_reads" %in% colnames(selected)) {
        readsPerUMI=mean(selected$num_reads/selected[[transcriptFeature]])
        strTitle=sprintf("%s, intronic>=%.2f\n%s Nuclei, %.1f reads/UMI, %d+ UMIs, medUMIs %s",
                         strTitlePrefix, minIntronic, comma_formatter(numSTAMPs), readsPerUMI, minNumUmis, comma_formatter(medianUMI))
    } else {

    }

    return (strTitle)

}

#' Plot the expression vs intronic content of the cells, colored by the fraction of UMIs removed by CellBender.
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
plotExpressionVsIntronic <- function(cell_features, densityCenters = NULL, intronicThreshold = NULL, title = "All cell barcodes", point_size=0.25, alpha=0.25, useCellBenderFeatures=TRUE) {
    if (!useCellBenderFeatures) {
        p=ggplot() + theme_void() #empty plot
        return (p)
    }

    #TO MAKE R CMD CHECK HAPPY
    x=y=num_transcripts=pct_intronic=frac_contamination=NULL;

    p <- ggplot(cell_features, aes(x = log10(num_transcripts), y = pct_intronic, color = frac_contamination)) +
        ggrastr::rasterize(geom_point(size = point_size, alpha=alpha),dpi=900) +
        scale_color_gradientn(colors = c("light blue", "blue", "red"), values = c(0,0.5,1)) +
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
            geom_point(data = data.frame(x = log10(densityCenters$empty[1]), y = densityCenters$empty[2]), aes(x = x, y = y), color = "red", size = 3) +
            geom_point(data = data.frame(x = log10(densityCenters$cell[1]), y = densityCenters$cell[2]), aes(x = x, y = y), color = "red", size = 3)
    }

    if (!is.null(intronicThreshold)) {
        p <- p +
            geom_hline(yintercept = intronicThreshold, linetype = "dotted", color = "red", linewidth = 1)
    }

    return(p)
}


#' Plot the expression vs intronic content of the cells
#' @param cell_features A data frame containing the cell features.
#' @param bounds_empty A list containing the bounds for the empty cells.
#' @param bounds_non_empty A list containing the bounds for the non-empty cells.
#' @param strTitleOverride The title of the plot.  Overrides a more generic title.
#' @param cex.axis The size of the axis text.
#' @param cex.lab The size of the axis labels.
#' @param cex.main The size of the main title.
#' @importFrom graphics par smoothScatter title rect
#' @noRd
plotCellTypeIntervals<-function (cell_features, bounds_empty, bounds_non_empty, strTitleOverride=NULL,
                                 cex.axis=0.6, cex.lab=0.7, cex.main=0.65) {

    par(mar=c(2,2,2,1), mgp=c(0.8,0.25,0),tck = -0.02)

    #how many cell barcodes in the training set?
    num_training_cells=length(which(cell_features$training_label_is_cell))

    smoothScatter(
        log10(cell_features$num_transcripts),
        cell_features$pct_intronic,
        xlim=log10_UMI_AXIS_RANGE_NEW,
        ylim=c(0,1),
        xlab="log10( UMIs )", ylab="intronic",
        main="", axes=TRUE, cex.axis=cex.axis, cex.lab=cex.lab
    )

    strTitle=paste("SVM Training -- initial nuclei exemplars selected [", num_training_cells, "]", sep="")
    #optionally override title.
    if (!is.null(strTitleOverride))
        strTitle=strTitleOverride

    title(main=strTitle, line=1, col.main="black", cex.main=cex.main)

    graphics::rect(bounds_empty$umi_lower_bound, bounds_empty$intronic_lower_bound, bounds_empty$umi_upper_bound, bounds_empty$intronic_upper_bound, border="red", lwd=3, lty=1)
    graphics::rect(bounds_non_empty$umi_lower_bound, bounds_non_empty$intronic_lower_bound, bounds_non_empty$umi_upper_bound, bounds_non_empty$intronic_upper_bound, border="green", lwd=3, lty=1)

    #reset par to defaults.  trying to save the original PAR and reset it was breaking.
    par(mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0), tck = NA)
}

plotSelectedCells<-function (cell_features_result, size = 0.25, alpha=0.25) {

    strTitle="Selected Nuclei"

    umi_min_threshold=min (log10(cell_features_result[cell_features_result$is_cell==T,][["num_transcripts"]]))
    intronic_min_threshold=min (cell_features_result[cell_features_result$is_cell==T,]$pct_intronic)

    #TO MAKE R CMD CHECK HAPPY
    num_transcripts=pct_intronic=is_cell=NULL;

    p <- ggplot(cell_features_result, aes(x = log10(num_transcripts), y = pct_intronic, color = is_cell)) +
        ggrastr::rasterize(geom_point(size = size, alpha=alpha),dpi=900) +
        labs(
            x = "log10(UMI)",
            y = "% Intronic",
            color = "Selected Nuclei"
        ) +
        ggtitle(strTitle) +
        scale_color_manual(values = c("TRUE" = "green", "FALSE" = "lightblue")) +
        coord_cartesian(xlim = log10_UMI_AXIS_RANGE_NEW) +
        theme_minimal() +
        guides(color = guide_legend(override.aes = list(size = 4, alpha=1)))

    return (p)
}

#changePar makes this plot compatable with the standard analysis cell selection pipeline that uses this plot in-line with standard graphics.

#' Plot Selected Cells Smooth Scatter
#'
#' This function creates a smooth scatter plot for selected cells, showing the relationship
#' between transcript counts and percent intronic features.
#'
#' @param cell_features A dataframe containing cell metadata with at least the columns
#'   specified by `transcriptFeature` and `pct_intronic`.
#' @param strTitlePrefix A prefix for the plot title.
#' @param transcriptFeature The column name for transcript counts (default: "num_transcripts").
#' @param changePar Logical; if TRUE, adjusts graphical parameters for the plot.
#' @param useCellBenderFeatures Logical; if FALSE, returns an empty ggplot object.
#'
#' @return A plot object or NULL.
#' @import ggplot2
#' @importFrom graphics par smoothScatter title abline
#' @noRd
plotSelectedCellsSmoothScatter <- function(cell_features, strTitlePrefix = "",
                                           transcriptFeature = "num_transcripts",
                                           changePar = TRUE, useCellBenderFeatures = TRUE) {
    if (!useCellBenderFeatures) {
        p <- ggplot2::ggplot() + ggplot2::theme_void() # empty plot
        return(p)
    }

    xlab <- "log10 ( UMIs )"
    strPrefix <- ""

    if (transcriptFeature == "num_retained_transcripts") {
        xlab <- "log10 ( UMIs post CBRB )"
    }

    cell_features_filtered <- cell_features[!is.na(cell_features$is_cell), ]
    #for cases where plotting remove background processed data.
    #number of transcripts should never be 0.
    cell_features_filtered <- cell_features_filtered[cell_features_filtered[[transcriptFeature]]>0,]

    strTitle <- getCellSelectionPlotTitle(cell_features_filtered,
                                          strTitlePrefix = strTitlePrefix,
                                          transcriptFeature = transcriptFeature)

    umi_min_threshold <- min(cell_features_filtered[cell_features_filtered$is_cell == TRUE, ][[transcriptFeature]])
    intronic_min_threshold <- min(cell_features_filtered[cell_features_filtered$is_cell == TRUE, ]$pct_intronic)

    if (changePar) {
        opar <- graphics::par(no.readonly = TRUE)
        graphics::par(mar = c(2, 2, 2, 1), mgp = c(0.8, 0.25, 0), tck = -0.02)
    }

    graphics::smoothScatter(
        log10(cell_features[[transcriptFeature]]),
        cell_features$pct_intronic,
        xlim = log10_UMI_AXIS_RANGE_NEW,
        ylim = c(0, 1),
        xlab = xlab,
        ylab = "intronic",
        axes = TRUE,
        cex.axis = 0.6,
        cex.lab = 0.7
    )

    graphics::title(main = strTitle, line = 0.25, col.main = "black", cex.main = 0.65)
    graphics::abline(v = log10(umi_min_threshold), col = "red", lty = 2)
    graphics::abline(h = intronic_min_threshold, col = "red", lty = 2)

    if (changePar) {
        graphics::par(opar)
    }
}

plotScaledTrainingDataFeatures<-function(trainingData) {
    # This is to avoid using dplyr and adding more dependencies.
    # Find the index of the label column
    label_col_index <- which(names(trainingData) == "training_label_is_cell")

    # Exclude the label column when generating feature names and values
    feature_columns <- trainingData[, -label_col_index]
    label_column <- trainingData[, label_col_index]

    # Create the long-format data frame
    longData <- data.frame(
        feature = rep(names(feature_columns), each = nrow(trainingData)),
        value = unlist(feature_columns),
        is_nuclei = rep(label_column, times = ncol(feature_columns))
    )

    # Plot
    #TO MAKE R CMD CHECK HAPPY
    feature=value=is_nuclei=NULL;

    p=ggplot(longData, aes(x = feature, y = value, fill = is_nuclei)) +
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
    return (p)
}

# plotCellProbability<- function(cell_features_result, title = "Assignment ") {
#     #TO MAKE R CMD CHECK HAPPY
#     num_transcripts=pct_intronic=is_cell_prob=NULL;
#
#     p <- ggplot(cell_features_result, aes(x = log10(num_transcripts), y = pct_intronic, color = is_cell_prob)) +
#         ggrastr::rasterize(geom_point(size = 0.25, alpha=0.25),dpi=900) +
#         labs(
#             x = "log10(UMI)",
#             y = "% Intronic",
#             color = "Cell Probability"
#         ) +
#         coord_cartesian(xlim = log10_UMI_AXIS_RANGE_NEW) +
#         ggtitle(title) +
#         theme_minimal()
#
#     return(p)
# }

plotCellProbabilities<-function (cell_features, strTitle="Nuclei Probability") {

    df=cell_features[which(cell_features$is_cell_prob>=0.5),]

    breaks <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
    colors <- c('red', 'purple', 'orange', 'blue', 'green')

    # Assign bins to is_cell_prob
    df$is_cell_prob_bin <- cut(df$is_cell_prob, breaks = breaks, include.lowest = TRUE)

    # Plot

    #TO MAKE R CMD CHECK HAPPY
    num_transcripts=pct_intronic=is_cell_prob_bin=NULL;

    p <- ggplot(df, aes(x = log10(num_transcripts), y = pct_intronic, color = is_cell_prob_bin)) +
        ggrastr::rasterize(geom_point(size = 0.5, alpha=0.25), dpi=900) +
        labs(
            x = "log10(UMI)",
            y = "% Intronic",
            color = "Nuclei Probability"
        ) +
        #coord_cartesian(xlim = log10_UMI_AXIS_RANGE_NEW) +
        ggtitle(strTitle) +
        theme_minimal() +
        scale_color_manual(values = colors) +
        guides(color = guide_legend(override.aes = list(size = 5, alpha=1)))

    return(p)
}



#######################
# GENE MODULE PLOTS
#######################

scatterPlotModuleScore <- function(cell_features, moduleName = "nuclei_gene_module_score", strTitle="", point_size = 0.25, alpha = 0.25) {
    #TO MAKE R CMD CHECK HAPPY
    num_transcripts=pct_intronic=NULL;

    if (!moduleName %in% colnames(cell_features)){
        stop(paste("Module name ", moduleName, " not found in cell features"))
    }

    p <- ggplot(cell_features, aes(x = log10(num_transcripts), y = pct_intronic, color = get(moduleName))) +
        ggrastr::rasterize(geom_point(size = point_size, alpha = alpha), dpi = 900) +
        scale_color_gradient2(low = "red", mid = "lightgrey", high = "blue", midpoint = 0) +
        labs(
            x = "log10(UMI)",
            y = "% Intronic",
            color = ""
        ) +
        coord_cartesian(xlim = log10_UMI_AXIS_RANGE_NEW) +
        ggtitle(moduleName) +
        theme_minimal() # +
    #custom_theme(title_size = 8, axis_title_size = 8, axis_text_size = 8)
    return(p)
}

plotModuleVsFracContam <- function(cell_features, moduleName = "nuclei_gene_module_score", point_size = 0.5, alpha = 0.5) {
    #TO MAKE R CMD CHECK HAPPY
    frac_contamination=pct_intronic=NULL;

    if (!moduleName %in% colnames(cell_features)){
        stop(paste("Module name ", moduleName, " not found in cell features"))
    }

    p <- ggplot(cell_features, aes(x = frac_contamination, y=get(moduleName), color = pct_intronic)) +
        ggrastr::rasterize(geom_point(size = point_size, alpha = alpha),dpi=900) +
        labs(
            x = "Fraction Contamination",
            y = moduleName,
            color = "% Intronic"
        ) +
        theme_minimal()
    return(p)
}

plotCellProbabilityConditionalCbrb<-function (cell_features_result, frac_contamination_threshold=1, useCellBenderFeatures=TRUE) {


    #if not using celbender, short circuit with an empty plot
    if (!useCellBenderFeatures) {
        return (ggplot() + theme_minimal())
    }

    #df=cell_features_result[cell_features_result$frac_contamination>=frac_contamination_threshold & cell_features_result$is_cell==T,]

    df <- cell_features_result[cell_features_result$frac_contamination >= frac_contamination_threshold,]

    # Filter out NAs
    df <- df[!is.na(df$is_cell_prob), ]

    # Calculate the number of cell barcodes
    numCellBarcodes <- dim(df[df$is_cell == TRUE,])[1]

    r=graphics::hist(df$is_cell_prob, breaks=seq(0, 1, by = 0.01), plot=F)

    dd=data.frame(countLog10=log10(r$counts+1), midpoint=r$mids)

    #TO MAKE R CMD CHECK HAPPY
    midpoint=countLog10=NULL;

    p <- ggplot(dd, aes(x = midpoint, y = countLog10)) +
        geom_bar(stat = "identity", width = 0.01) +
        geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", size = 1) +
        labs(x = "Cell Probability", y = "Log10(Count + 1)") +
        ggtitle(paste("SVM Probability for CBRB empty cell barcodes\n[", numCellBarcodes, "] barcode selected as nuclei")) +
        theme_minimal() +
        annotation_logticks(sides = "l")

    return (p)

}






##########################################
# OTHER CODE - to be organized or removed.
###########################################

