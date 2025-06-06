####################
# HIGH LEVEL BOUNDS FUNCTIONS
#####################
#' Find exemplars for nuclei and empty droplets
#'
#' @param cell_features The cell features data frame
#' @param max_umis_empty Exclude cell barcodes from analysis with fewer than
#'   these many UMIs.  This threshold should exclude noisy barcodes, but not
#'   exclude the empty droplet cloud.
#' @param useCellBenderFeatures When true, the CellBender remove-background
#'   feature frac_contamination is used for cell selection.
#' @param forceTwoClusterSolution When true, the function will attempt to find a
#'   solution with two clusters.  This may be useful when the data is overloaded
#'   and breaks the normal assumptions, but may find suboptimal solutions for
#'   other data sets.
#' @param verbose Print verbose output to log
#' @return A list containing the bounds for the empty droplet and nuclei
#'   exemplars
#' @import logger
#' @noRd
findTrainingDataBounds <- function(cell_features, max_umis_empty = 50,
    useCellBenderFeatures = TRUE, forceTwoClusterSolution = FALSE,
    verbose = FALSE) {

    # If using CellBender features, use the specific CBRB method.
    if (useCellBenderFeatures) {
        return(findTrainingDataBoundsCBRB(cell_features,
            max_umis_empty = max_umis_empty))

    }

    # A lower UMI threshold is used for some cases.
    max_umis_empty_off <- 20

    # If forcing a two-cluster solution, try to find separation between the two
    # highest peaks.
    if (forceTwoClusterSolution) {
        return(findTwoClusterSolution(cell_features, max_umis_empty,
            max_umis_empty_off, verbose))
    } else {
        # Otherwise, use a more general approach with empty droplet fraction
        # constraints.
        return(findDefaultSolution(cell_features, max_umis_empty,
            max_umis_empty_off, verbose))
    }
}

findTwoClusterSolution <- function(cell_features, max_umis_empty,
    max_umis_empty_off, verbose) {
    # With the default UMI filter, look for the separation between the two
    # highest peaks.
    logSelectionProcess("Separate Two Highest Peaks", max_umis_empty, verbose)
    defaultPitBetween <- findTrainingDataBoundsDefaultIterative(
        cell_features, max_umis_empty = max_umis_empty,
        method = "PitBetweenHighestPeaks", verbose = verbose)

    # Without a UMI filter, this method works better in high-noise experiments.
    logSelectionProcess("Separate Two Highest Peaks No Filter",
    max_umis_empty_off, verbose)

    noUmiFilterPitBetween <- findTrainingDataBoundsDefaultIterative(
        cell_features, max_umis_empty = max_umis_empty_off,
        method = "PitBetweenHighestPeaks", verbose = verbose)

    results <- list(`Separate Two Highest Peaks` = defaultPitBetween,
        `Separate Two Highest Peaks No Filter` = noUmiFilterPitBetween)

    # In forced two-cluster solutions, do not enforce a minimum fraction of
    # empty droplets.
    return(selectBestTrainingBoundsModel(results))
}

findDefaultSolution <- function(cell_features, max_umis_empty,
    max_umis_empty_off, verbose) {

    # The default approach finds the pit after the highest peak.
    logSelectionProcess("Pit After Highest Peak", max_umis_empty, verbose)
    default <- findTrainingDataBoundsDefaultIterative(
        cell_features, max_umis_empty = max_umis_empty, verbose = verbose)

    # Without a UMI filter, this is useful for low-noise experiments.
    logSelectionProcess("Pit After Highest Peak No Filter",
        max_umis_empty_off, verbose)

    noUmiFilter <- findTrainingDataBoundsDefaultIterative(
        cell_features, max_umis_empty = max_umis_empty_off, verbose = verbose)

    results <- list(`Default Selection` = default,
        `No UMI Filter Selection` = noUmiFilter)

    # If not forcing a two-cluster model, enforce a minimum fraction of empty
    # droplets!
    return(selectBestTrainingBoundsModel(results,
        minEmptyDropletFraction = 0.2))
}

logSelectionProcess <- function(method, umi_filter, verbose) {
    # Logs the current method and the UMI threshold being used.
    msg <- paste("Exemplar selection [",
        method, "] with UMI filter [", umi_filter, "]")

    log_info(msg)
}

#' Select the best model from a list
#'
#' This function selects the best model from a list of results.  If the minimum
#' empty droplet fraction is set, it will only consider models that have a
#' minimum fraction of empty droplets. If no models meet the minimum empty
#' droplet fraction, it will select the best model from the original list.
#'
#' @param results A list of results from findTrainingDataBoundsDefaultIterative
#' @param minEmptyDropletFraction The minimum fraction of empty droplets that
#'   must be present in the training data.
#' @return The best model from the list
#' @import logger
#' @noRd
selectBestTrainingBoundsModel <- function(results,
    minEmptyDropletFraction = NULL) {

    # To make bioconductor happy use vapply.
    silhouettes <-
        vapply(results, function(res) res$best_silhouette, numeric(1))
    numEmpty <- vapply(results, function(res) res$numEmpty, numeric(1))
    numNonEmpty <- vapply(results, function(res) res$numNonEmpty, numeric(1))

    # the row names of the data frame are the model names.
    df <- data.frame(silhouettes = silhouettes,
        numEmpty = numEmpty,
        numNonEmpty = numNonEmpty)

    # remove NA silhouette results.  this explicitly ignores NAs, in
    # cases where the model failed to find good bounds.
    df <- df[!is.na(df$silhouettes), ]
    best_name <- rownames(df[which.max(df$silhouettes), ])

    # if maxEmptyDropletFraction is not null, drop any results that
    # have a low number of empty droplets.
    if (!is.null(minEmptyDropletFraction)) {
        df <- df[df$numEmpty/df$numNonEmpty > minEmptyDropletFraction,]
    }

    # if there are results left with reasonable numbers of empty/non
    # empties, keep the best silhouette score.  if there are no
    # results left, keep the best silhouette score from the original
    # list.
    if (nrow(df) > 0) {
        best_name <- rownames(df[which.max(df$silhouettes), ])
    }

    # Log the selected result
    msg<-paste("Selected result [",best_name,"] with best silhouette score [",
        round(max(silhouettes), 3), "]", sep = "")
    log_info(msg)
    # label the results that have the best initialization
    result <- results[[best_name]]
    result$method <- best_name
    return(result)
}


#############################
# SIMPLE CBRB
#############################
findTrainingDataBoundsCBRB <- function(cell_features, max_umis_empty = 50,
    maxContaminationThreshold = 0.1) {
    # the interval for empty cells for training data explicitly avoid
    # searching in the very low UMI area, which are likely tons of
    # PCR error barcodes.
    cell_features_empty <- cell_features[cell_features$frac_contamination == 1 &
        cell_features$num_transcripts > max_umis_empty, ]

    bounds_empty <- getHighestDensityIntervalsEnforcedSmoothing(
        cell_features_empty, pctDensity = 75)

    # the intervals for nuclei examples
    b <- selectNucleiExemplarBounds(cell_features, maxContaminationThreshold,
        max_umis_empty = max_umis_empty, initialDensity = 95,
        bounds_empty = bounds_empty, extendCellSelectionBounds = TRUE)

    bounds_non_empty <- b$bounds_non_empty
    bounds_non_empty_extended <- b$bounds_non_empty_extended

    result <- list(bounds_empty = bounds_empty,
        bounds_non_empty = bounds_non_empty_extended)

    return(result)
}



################################
# BOUNDS DISCOVERY WTIHOUT CBRB
#################################

#' Find the training data bounds using an iterative approach
#'
#' This replaces findTrainingDataBoundsDefault with an iterative approach to
#' find the training data bounds. This function performs an iterative process to
#' select the initial empty/nuclei clusters by adjusting the density smoothing
#' factor, computing UMI thresholds, and optimizing based on silhouette scores.
#' The function applies a grid search over the density adjustments to identify
#' the optimal UMI threshold for separating nuclei and empty droplets.
#'
#' @param cell_features A dataframe containing an entry for each cell barcode
#'   with columns num_transcripts and pct_intronic.
#' @param max_umis_empty Exclude cell barcodes from analysis with fewer than
#'   this many UMIs. This threshold should exclude noisy barcodes, but not
#'   exclude the empty droplet cloud.
#' @param method A character string specifying the method to use for selecting
#'   the UMI threshold. Options are: PitAfterHighestPeak (default) and
#'   PitBetweenHighestPeaks. The former selects the UMI threshold after the
#'   highest peak in the density plot, while the latter selects the UMI
#'   threshold between the two highest peaks in the density plot.
#' @param smoothingMultiple A numeric multiplier that controls the rate of
#'   smoothing adjustment. Default is 1.1.
#' @param max_iterations Maximum number of iterations for adjusting the
#'   smoothing factor. Default is 10.
#' @param early_termination_threshold Numeric value representing the percentage
#'   drop in silhouette score for early termination. Default is 0.05.
#' @param verbose Print verbose output to log
#' @return A list containing A dataframe with columns smoothingMultipliers,
#'   umiThresholds, and silhouetteScores for each iteration. The bounds found
#'   for the maximum silhouette score and the maximum silhouette score.
#'
#' @import logger
#' @noRd
findTrainingDataBoundsDefaultIterative <- function(cell_features,
    max_umis_empty = 50,
    method = c("PitAfterHighestPeak", "PitBetweenHighestPeaks"),
    smoothingMultiple = 1.1, max_iterations = 10,
    early_termination_threshold = 0.05, verbose = TRUE) {
    method <- match.arg(method)
    log_info("Starting iterative smoothing with method [", method, "]")
    df_filtered <-
        cell_features[cell_features$num_transcripts >= max_umis_empty, ]
    x <- log10(df_filtered$num_transcripts + 1)
    best_silhouette <- -Inf; smoothingMultiplier <- 1
    best_bounds <- best_umi_threshold <- NULL
    results <- data.frame(smoothingMultiplier = numeric(0),
        umiThreshold = numeric(0), silhouetteScore = numeric(0))
    for (iter in seq_len(max_iterations)) {
        denRange <- seq(0.25, 3, 0.25)^smoothingMultiplier
        umiThreshold <- determineUMIThreshold(x, method, denRange)
        if (is.na(umiThreshold)) {
            results <- rbind(results, c(smoothingMultiplier, NA, NA))
            next
        }
        bounds <- findTrainingDataBoundsDefault(df_filtered, max_umis_empty,
            umiThresholdOverride = umiThreshold, verbose = FALSE)
        if (is.null(bounds)) {
            results <- rbind(results, c(smoothingMultiplier, umiThreshold, NA))
            next
        }
        cell_features_labeled <- labelTrainingData(df_filtered,
            bounds$bounds_empty,  bounds$bounds_non_empty, NULL,
            useCellBenderFeatures = FALSE, verbose = FALSE)
        silhouetteResult <- calculate_silhouette(cell_features_labeled,
            showPlot = FALSE, verbose = FALSE)
        score <- silhouetteResult$mean_silhouette
        results <- rbind(results, c(smoothingMultiplier, umiThreshold, score))
        if (!is.na(score) && score > best_silhouette) {
            best_silhouette <- score
            best_bounds <- bounds
            best_umi_threshold <- umiThreshold
        }
        threshold<-(1 - early_termination_threshold) * best_silhouette
        if (!is.na(score) && score < threshold) {
            msg<-"Early termination: silhouette score dropped significantly"
            log_info(msg)
            break
        }
        smoothingMultiplier <- smoothingMultiplier * smoothingMultiple
    }
    return(finalizeTrainingResults(df_filtered, best_silhouette,
        best_umi_threshold, best_bounds, results))
}

determineUMIThreshold <- function(x, method, denRange) {
    if (method == "PitAfterHighestPeak") {
        return(PitAfterHighestPeakWithGridSearch(x, adjust = denRange))
    }
    if (method == "PitBetweenHighestPeaks") {
        return(PitBetweenHighestPeaksWithGridSearch(x, adjust = denRange))
    }
    stop("Invalid method specified.")
}

finalizeTrainingResults <- function(df_filtered, best_silhouette,
    best_umi_threshold, best_bounds, results) {
    if (is.null(best_umi_threshold)) {
        log_info("Unable to find good initialization. Returning Empty Bounds")
        return(list(best_silhouette = NA, best_umi_threshold = NA,
            bounds_empty = NA, bounds_non_empty = NA, resultDF = results,
            numEmpty = NA, numNonEmpty = NA))
    }

    cell_features_labeled <-
        labelTrainingData(df_filtered, best_bounds$bounds_empty,
        best_bounds$bounds_non_empty, NULL,
        useCellBenderFeatures = FALSE, verbose = FALSE)

    numEmpty <- sum(!is.na(cell_features_labeled$training_label_is_cell) &
        !cell_features_labeled$training_label_is_cell)
    numNonEmpty <- sum(!is.na(cell_features_labeled$training_label_is_cell) &
        cell_features_labeled$training_label_is_cell)

    return(list(best_silhouette = best_silhouette,
        best_umi_threshold = best_umi_threshold,
        bounds_empty = best_bounds$bounds_empty,
        bounds_non_empty = best_bounds$bounds_non_empty,
        resultDF = results, numEmpty = numEmpty, numNonEmpty = numNonEmpty))
}

############
# FIND DEFAULT TRAINING DATA BOUNDS
###########

#' Find the training data bounds using only the \% intronic and num_transcripts
#' features
#'
#' This function uses a turning point method to find the largest peak and a peak
#' after that to define a linear separation of the data into the set that likely
#' contains empty droplets and the part of the data that contains non-empty
#' droplets. It then selects areas of maximum density within each of the two
#' partitions.
#'
#' @param cell_features A dataframe containing the dataset with at least columns
#'   num_transcripts and pct_intronic.
#' @param max_umis_empty Exclude cell barcodes from analysis with fewer than
#'   this many UMIs. This is used to exclude very small barcodes that are noise
#'   and may confuse the initial turning point threshold.
#' @param umiThresholdOverride A numeric value to override the UMI threshold. If
#'   NULL, the function will find the UMI threshold using the turning point
#'   method. If not NULL, use an initial estimate provided by another means.
#' @param verbose Logical. if TRUE, prints verbose output to the log.
#'
#' @return A list with the training data bounds.
#' @noRd
findTrainingDataBoundsDefault <- function(cell_features, max_umis_empty = 50,
    umiThresholdOverride = NULL, verbose = TRUE) {

    df <- cell_features[cell_features$num_transcripts >= max_umis_empty, ]

    umiThreshold <- if (!is.null(umiThresholdOverride)) umiThresholdOverride
    else PitAfterHighestPeakWithGridSearch(log10(df$num_transcripts))

    if (verbose)
        log_info("Initial UMI Threshold [", round(umiThreshold, 3), "]")

    df_empty <- df[log10(df$num_transcripts) < umiThreshold, ]
    bounds_empty <- getEmptyCellsByDensity(df_empty, "pct_intronic",
        pctDensity = 75, showPlot = FALSE, verbose = verbose)

    early_exit <- checkEarlyExit(df_empty, bounds_empty, umiThreshold, verbose)
    if (!is.null(early_exit)) return(early_exit)

    df_non_empty <- df[log10(df$num_transcripts) > umiThreshold, ]
    early_exit <- checkEarlyExit(df_non_empty, bounds_empty,
        umiThreshold, verbose)
    if (!is.null(early_exit)) return(early_exit)

    bounds_non_empty_intronic <- getHighestDensityIntervalsEnforcedSmoothing(
        df_non_empty, yAxisFeature = "pct_intronic", pctDensity = 75,
        maxPeaksExpected = 1, showPlot = FALSE
    )
    df_filtered <- df_non_empty[
        df_non_empty$pct_intronic >=
            bounds_non_empty_intronic$intronic_lower_bound &
        df_non_empty$pct_intronic <=
            bounds_non_empty_intronic$intronic_upper_bound,]

    bounds_non_empty_transcripts <- getHighestDensityIntervalsEnforcedSmoothing(
        df_filtered, yAxisFeature = "num_transcripts",
        pctDensity = 75, showPlot = FALSE
    )
    bounds_non_empty_transcripts$umi_upper_bound <-
        as.numeric(stats::quantile(log10(df_filtered$num_transcripts), 0.995))

    bounds_non_empty <- data.frame(
        umi_lower_bound = bounds_non_empty_transcripts$umi_lower_bound,
        umi_upper_bound = bounds_non_empty_transcripts$umi_upper_bound,
        intronic_lower_bound = bounds_non_empty_intronic$intronic_lower_bound,
        intronic_upper_bound = bounds_non_empty_intronic$intronic_upper_bound
    )

    return(list(bounds_empty = bounds_empty,
        bounds_non_empty = bounds_non_empty))
}

checkEarlyExit <- function(df, bounds_empty, umiThreshold, verbose) {
    if (nrow(df) == 0 || is.null(bounds_empty)) {
        if (verbose) log_warn("No data left after selecting empty barcodes.")
        bounds_empty$umi_upper_bound <- umiThreshold
        return(list(
            bounds_empty = bounds_empty,
            bounds_non_empty = data.frame(
                umi_lower_bound = NA, umi_upper_bound = NA,
                intronic_lower_bound = NA, intronic_upper_bound = NA
            )
        ))
    }
    return(NULL)
}

# if there are multiple peaks, only keep the largest peak.
getEmptyCellsByDensity <- function(cell_features, yAxisFeature = "pct_intronic",
    pctDensity = 75, showPlot = FALSE, verbose = FALSE) {
    probList <- unique(c(25, 50, 75, 90, 95, 99, pctDensity))
    probList <- probList[probList <= pctDensity]

    # Do this by density.  This distribution is bimodal (empty/nuclei.)
    x <- log10(cell_features$num_transcripts + 1)
    zHDR <- getBoundsByDensity(x, probList = probList, pctDensity = pctDensity,
        maxPeaksExpected = 2, showPlot = showPlot)

    # which interval contains the majority of the data?
    zHDR_intervals <- matrix(zHDR$intervals, ncol = 2, byrow = TRUE)

    # copied from the internals of hdrcde::hdr.den we use the same
    # bandwidth as used in getBoundsByDensity.
    d <- density(x, bw = zHDR$bw, n = 1001)
    maxDensityX <- d$x[which.max(d$y)]

    # Step 2: Calculate the fraction of scores in X within each
    # interval
    containsInterval <- function(interval, maxDensityX) {
        return(maxDensityX >= interval[1] & maxDensityX <= interval[2])
    }

    # Step 3: Identify the interval with the highest density
    max_idx <- which(apply(zHDR_intervals, 1, containsInterval, maxDensityX))
    # this step can fail when the initial cutpoint is poor.
    if (length(max_idx) == 0) {
        if (verbose) {
            msg<-paste("Unable to select the empty cell barcode highest",
                "density interval.")
            log_warn(msg)
        }
        return(NULL)
    }
    umiBounds <- zHDR_intervals[max_idx, ]

    # select data from this interval and find the %intronic.
    df <- cell_features[log10(cell_features$num_transcripts) >= umiBounds[1] &
        log10(cell_features$num_transcripts) <= umiBounds[2], ]
    pctDensity <- 95
    probList <- unique(c(25, 50, 75, 90, 95, 99, pctDensity))
    probList <- probList[probList <= pctDensity]
    yBounds <- getBoundsByDensity(df[[yAxisFeature]], probList, pctDensity,
        maxPeaksExpected = 1, showPlot = showPlot)$intervals
    result <- data.frame(umi_lower_bound = umiBounds[1],
        umi_upper_bound = umiBounds[2], intronic_lower_bound = yBounds[1],
        intronic_upper_bound = yBounds[2])
    return(result)
}
