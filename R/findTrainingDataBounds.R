####################
# HIGH LEVEL BOUNDS FUNCTIONS
#####################

#' Find exemplar bounds for SVM training
#'
#' This function selects exemplar barcode regions for downstream SVM training.
#' The initialization strategy is controlled by `useCBRBInitialization`.
#'
#' When `useCBRBInitialization` is `TRUE`, exemplar bounds are selected with
#' the CellBender remove-background initialization path. This path is intended
#' for empty-vs-nucleus initialization and does not define debris exemplars.
#'
#' When `useCBRBInitialization` is `FALSE`, exemplar bounds are selected with
#' the default density-based initialization. The algorithm first estimates a UMI
#' threshold that separates likely-empty and likely-nucleus barcode partitions.
#' It then selects high-density exemplar regions within each partition using
#' log10 UMI counts and pct_intronic. This path can also define debris
#' exemplars.
#'
#' In the default density-based path, low-intronic debris can fall in the
#' high-UMI partition and distort nucleus exemplar selection. To reduce this
#' failure mode, the algorithm applies an intronic floor before estimating the
#' nucleus exemplar bounds. The observed floor is computed from the empty
#' droplet exemplar bounds:
#'
#'   observed_intronic_floor =
#'       intronic_floor_fraction * bounds_empty$intronic_lower_bound
#'
#' This observed floor is then capped by `debris_pct_intronic_prior`:
#'
#'   applied_intronic_floor =
#'       min(observed_intronic_floor, debris_pct_intronic_prior)
#'
#' The high-UMI candidate nucleus partition is filtered to barcodes with
#' pct_intronic greater than or equal to `applied_intronic_floor` before nucleus
#' bounds are estimated. This allows the empty droplet distribution to define a
#' data-driven lower bound when the empty droplet intronic floor is low, while
#' preventing high-intronic empty droplet clusters from forcing an overly strict
#' nucleus intronic cutoff.
#'
#' @param cell_features A data frame containing barcode-level features.
#' @param max_umis_empty Exclude barcodes with fewer than this many UMIs from
#'   exemplar-bound detection. This threshold should remove noisy low-count
#'   barcodes while retaining the empty droplet cloud.
#' @param useCBRBInitialization Logical scalar. If `TRUE`, use the CBRB
#'   initialization path to estimate empty and nucleus exemplar bounds. If
#'   `FALSE`, use the default density-based initialization, which can also
#'   identify debris exemplars.
#' @param forceTwoClusterSolution Logical scalar. If `TRUE`, attempt to find a
#'   two-cluster solution by separating the two highest-density peaks. This can
#'   help overloaded datasets but may be suboptimal for typical datasets.
#' @param intronic_floor_fraction Numeric scalar. Multiplier applied to the
#'   empty droplet intronic lower bound when computing the observed intronic
#'   floor used to remove low-intronic debris from the candidate nucleus
#'   partition.
#' @param debris_pct_intronic_prior Numeric scalar. Prior upper bound for the
#'   pct_intronic range of debris-like barcodes. The applied intronic floor is
#'   the smaller of this value and the adjusted empty-derived intronic lower
#'   bound.
#' @param use2DTrainingRefinement Logical scalar. If `TRUE`, refine default
#'   density-based empty and nucleus exemplar selections with a two-dimensional
#'   HDR component selection. This is experimental and is not applied to debris.
#' @param verbose Logical scalar. If `TRUE`, emit diagnostic log messages.
#'
#' @return A list containing exemplar bounds and selected training barcode
#'   vectors. The default density-based path may include `bounds_debris` and
#'   `training_debris_barcodes`; the CBRB initialization path does not define
#'   debris exemplars.
#'
#' @import logger
#' @noRd
#'
findTrainingDataBounds <- function(
  cell_features, max_umis_empty = 50,
  useCBRBInitialization = TRUE, forceTwoClusterSolution = FALSE,
  intronic_floor_fraction = 0.9, debris_pct_intronic_prior = 0.25,
  use2DTrainingRefinement = FALSE, verbose = FALSE
) {
  # If using CBRB initialization, use the CBRB-specific bounds method.
  if (useCBRBInitialization) {
    return(findTrainingDataBoundsCBRB(cell_features,
      max_umis_empty = max_umis_empty
    ))
  }

  # A lower UMI threshold is used for some cases.
  max_umis_empty_off <- 20

  # If forcing a two-cluster solution, try to find separation between the two
  # highest peaks.
  if (forceTwoClusterSolution) {
    return(findTwoClusterSolution(
      cell_features, max_umis_empty,
      max_umis_empty_off, intronic_floor_fraction,
      debris_pct_intronic_prior, use2DTrainingRefinement, verbose
    ))
  } else {
    # Otherwise, use a more general approach with empty droplet fraction
    # constraints.
    return(findDefaultSolution(
      cell_features, max_umis_empty,
      max_umis_empty_off, intronic_floor_fraction,
      debris_pct_intronic_prior, use2DTrainingRefinement, verbose
    ))
  }
}

findTwoClusterSolution <- function(
  cell_features, max_umis_empty,
  max_umis_empty_off, intronic_floor_fraction,
  debris_pct_intronic_prior, use2DTrainingRefinement, verbose
) {
  # With the default UMI filter, look for the separation between the two
  # highest peaks.
  logSelectionProcess("Separate Two Highest Peaks", max_umis_empty, verbose)
  defaultPitBetween <- findTrainingDataBoundsDefaultIterative(
    cell_features,
    max_umis_empty = max_umis_empty,
    method = "PitBetweenHighestPeaks",
    intronic_floor_fraction = intronic_floor_fraction,
    debris_pct_intronic_prior = debris_pct_intronic_prior,
    use2DTrainingRefinement = use2DTrainingRefinement,
    verbose = verbose
  )

  # Without a UMI filter, this method works better in high-noise experiments.
  logSelectionProcess(
    "Separate Two Highest Peaks No Filter",
    max_umis_empty_off, verbose
  )

  noUmiFilterPitBetween <- findTrainingDataBoundsDefaultIterative(
    cell_features,
    max_umis_empty = max_umis_empty_off,
    method = "PitBetweenHighestPeaks",
    intronic_floor_fraction = intronic_floor_fraction,
    debris_pct_intronic_prior = debris_pct_intronic_prior,
    use2DTrainingRefinement = use2DTrainingRefinement,
    verbose = verbose
  )

  results <- list(
    `Separate Two Highest Peaks` = defaultPitBetween,
    `Separate Two Highest Peaks No Filter` = noUmiFilterPitBetween
  )

  # In forced two-cluster solutions, do not enforce a minimum fraction of
  # empty droplets.
  return(selectBestTrainingBoundsModel(results))
}

findDefaultSolution <- function(
  cell_features, max_umis_empty,
  max_umis_empty_off, intronic_floor_fraction,
  debris_pct_intronic_prior, use2DTrainingRefinement, verbose
) {
  # The default approach finds the pit after the highest peak.
  logSelectionProcess("Pit After Highest Peak", max_umis_empty, verbose)
  default <- findTrainingDataBoundsDefaultIterative(
    cell_features,
    max_umis_empty = max_umis_empty,
    intronic_floor_fraction = intronic_floor_fraction,
    debris_pct_intronic_prior = debris_pct_intronic_prior,
    use2DTrainingRefinement = use2DTrainingRefinement,
    verbose = verbose
  )

  # Without a UMI filter, this is useful for low-noise experiments.
  logSelectionProcess(
    "Pit After Highest Peak No Filter",
    max_umis_empty_off, verbose
  )

  noUmiFilter <- findTrainingDataBoundsDefaultIterative(
    cell_features,
    max_umis_empty = max_umis_empty_off,
    intronic_floor_fraction = intronic_floor_fraction,
    debris_pct_intronic_prior = debris_pct_intronic_prior,
    use2DTrainingRefinement = use2DTrainingRefinement,
    verbose = verbose
  )

  results <- list(
    `Default Selection` = default,
    `No UMI Filter Selection` = noUmiFilter
  )

  # If not forcing a two-cluster model, enforce a minimum fraction of empty
  # droplets!
  return(selectBestTrainingBoundsModel(results,
    minEmptyDropletFraction = 0.2
  ))
}

logSelectionProcess <- function(method, umi_filter, verbose) {
  # Logs the current method and the UMI threshold being used.
  msg <- paste(
    "Exemplar selection [",
    method, "] with UMI filter [", umi_filter, "]"
  )

  log_info(msg)
}

#' Select the best model from a list
#'
#' This function selects the best candidate training-bound result. Results are
#' compared using the stability-regularized selection score. If
#' `minEmptyDropletFraction` is supplied, the function first tries to select
#' among results satisfying `numEmpty / numNonEmpty > minEmptyDropletFraction`.
#' If no result satisfies that filter, the function falls back to all candidate
#' results.
#'
#' @param results A named list of results from
#'     `findTrainingDataBoundsDefaultIterative()`.
#' @param minEmptyDropletFraction Numeric scalar or `NULL`. If supplied,
#'     candidates satisfying this empty-to-non-empty ratio are preferred.
#'
#' @return The selected training-bound result.
#' @noRd
selectBestTrainingBoundsModel <- function(
  results,
  minEmptyDropletFraction = NULL
) {
  makeResultSummaryRow <- function(resultName, result) {
    data.frame(
      resultName = resultName,
      selectionScore = result$best_selection_score,
      silhouette = result$best_silhouette,
      stability = result$best_stability,
      numEmpty = result$numEmpty,
      numNonEmpty = result$numNonEmpty,
      stringsAsFactors = FALSE
    )
  }

  resultDF <- do.call(
    rbind,
    Map(makeResultSummaryRow, names(results), results)
  )
  candidateDF <- resultDF

  if (!is.null(minEmptyDropletFraction)) {
    filteredDF <- resultDF[
      resultDF$numEmpty / resultDF$numNonEmpty >
        minEmptyDropletFraction,
    ]

    if (nrow(filteredDF) > 0) {
      candidateDF <- filteredDF
    }
  }

  bestRow <- candidateDF[which.max(candidateDF$selectionScore), ]
  bestName <- bestRow$resultName

  logger::log_info(
    paste0(
      "Selected result [", bestName, "] with selection score [",
      round(bestRow$selectionScore, 3), "] silhouette [",
      round(bestRow$silhouette, 3), "] stability [",
      round(bestRow$stability, 3), "]"
    )
  )

  result <- results[[bestName]]
  result$method <- bestName
  result
}
#############################
# SIMPLE CBRB
#############################
findTrainingDataBoundsCBRB <- function(
  cell_features, max_umis_empty = 50,
  maxContaminationThreshold = 0.1
) {
  # the interval for empty cells for training data explicitly avoid
  # searching in the very low UMI area, which are likely tons of
  # PCR error barcodes.
  cell_features_empty <- cell_features[cell_features$frac_contamination == 1 &
    cell_features$num_transcripts > max_umis_empty, ]

  bounds_empty <- getHighestDensityIntervalsEnforcedSmoothing(
    cell_features_empty,
    pctDensity = 75
  )

  # the intervals for nuclei examples
  b <- selectNucleiExemplarBounds(cell_features, maxContaminationThreshold,
    max_umis_empty = max_umis_empty, initialDensity = 95,
    bounds_empty = bounds_empty, extendCellSelectionBounds = TRUE
  )

  bounds_non_empty <- b$bounds_non_empty
  bounds_non_empty_extended <- b$bounds_non_empty_extended

  result <- list(
    bounds_empty = bounds_empty,
    bounds_non_empty = bounds_non_empty_extended
  )

  return(result)
}


################################
# BOUNDS DISCOVERY WITHOUT CBRB
#################################

#' Find the training data bounds using iterative smoothing
#'
#' This function searches over increasingly smoothed UMI-density estimates and
#' selects exemplar bounds for empty droplets and nuclei. Each smoothing value
#' generates candidate training-bound solutions that are scored by silhouette and
#' regularized by local stability across nearby smoothing values. The procedure
#' is designed to keep the UMI-threshold estimate robust to low-intronic debris
#' while preserving plausible empty-droplet and nucleus exemplar bounds.
#'
#' For each smoothing value, the algorithm runs the following steps:
#'
#' \enumerate{
#'   \item Estimate an initial UMI threshold from the full filtered UMI
#'     distribution. This threshold is provisional and is used only to obtain an
#'     initial empty-droplet estimate.
#'   \item Use the provisional threshold to estimate empty-droplet bounds with
#'     `findTrainingDataBoundsDefault()`.
#'   \item Use the provisional empty-droplet UMI upper bound to remove likely
#'     high-UMI debris from the vector used for UMI-threshold learning. A barcode
#'     is removed from threshold learning only if it is to the right of the
#'     provisional empty-droplet UMI bound and below `debris_pct_intronic_prior`.
#'     The full filtered data frame is still retained for final exemplar
#'     selection.
#'   \item Estimate a refined UMI threshold from the debris-filtered UMI vector.
#'   \item Recompute final empty and nucleus exemplar bounds from the full
#'     filtered data frame using the refined UMI threshold.
#'   \item Evaluate two candidate debris floors for the final nucleus partition:
#'     the default empty-derived floor and the fixed
#'     `debris_pct_intronic_prior` floor.
#'   \item Label training exemplars from each candidate solution and compute the
#'     silhouette score.
#' }
#'
#' After all smoothing values are evaluated, candidate solutions are assigned
#' local stability scores based on how consistently their empty and nucleus
#' exemplar sets are preserved across nearby smoothing values with the same
#' debris-floor source. The final candidate is selected using the
#' stability-regularized score.
#'
#' @param cell_features A data frame containing barcode-level features.
#' @param max_umis_empty Exclude cell barcodes from analysis with fewer than
#'     this many UMIs before estimating UMI thresholds.
#' @param method Method used to select the UMI threshold.
#' @param smoothingMultiple Multiplicative factor used to increase smoothing
#'     across iterations.
#' @param max_iterations Maximum number of smoothing iterations.
#' @param early_termination_threshold Fractional silhouette drop used for early
#'     stopping. If `NULL`, all smoothing candidates are evaluated.
#' @param stabilityLambda Maximum silhouette penalty applied to a completely
#'     unstable candidate.
#' @param stabilityNeighborWindow Number of nearby candidates on each side used
#'     to estimate local stability within each debris-floor source.
#' @param intronic_floor_fraction Fraction of the empty-droplet intronic lower
#'     bound used to compute the empty-derived floor for the high-UMI candidate
#'     nucleus partition.
#' @param debris_pct_intronic_prior Numeric scalar. Prior expectation for the
#'     upper pct_intronic range of debris-like barcodes. This value is used both
#'     when filtering high-UMI debris from UMI-threshold learning and as an
#'     alternative candidate floor for final nucleus exemplar selection.
#' @param verbose Logical scalar. If `TRUE`, log progress.
#'
#' @return A list containing selected bounds, selection metrics, diagnostic
#'     values for the selected candidate, and the per-iteration diagnostic data
#'     frame.
#' @noRd
findTrainingDataBoundsDefaultIterative <- function(
  cell_features,
  max_umis_empty = 50,
  method = c("PitAfterHighestPeak", "PitBetweenHighestPeaks"),
  smoothingMultiple = 1.1,
  max_iterations = 15,
  early_termination_threshold = NULL,
  stabilityLambda = 0.02,
  stabilityNeighborWindow = 2,
  intronic_floor_fraction = 1,
  debris_pct_intronic_prior = 0.25,
  use2DTrainingRefinement = FALSE,
  verbose = TRUE
) {
  method <- match.arg(method)

  log_info("Starting iterative smoothing with method [", method, "]")

  # Remove very small barcodes before threshold learning. These barcodes are
  # usually sequencing/PCR noise and can distort the first UMI density estimate.
  df_filtered <- cell_features[cell_features$num_transcripts >= max_umis_empty, ]
  x <- log10(df_filtered$num_transcripts + 1)

  best_silhouette <- -Inf
  smoothingMultiplier <- 1
  candidate_results <- list()
  resultDF <- makeEmptySmoothingResultDF()

  for (iter in seq_len(max_iterations)) {
    denRange <- seq(0.25, 3, 0.25)^smoothingMultiplier

    # First pass: estimate a rough UMI threshold from the full UMI density.
    # This pass is used only to obtain a provisional empty-droplet estimate.
    initialUmiThreshold <- determineUMIThreshold(x, method, denRange)

    if (is.na(initialUmiThreshold)) {
      resultDF <- addFailedSmoothingResult(
        resultDF = resultDF,
        smoothingMultiplier = smoothingMultiplier,
        umiThreshold = NA_real_
      )

      smoothingMultiplier <- smoothingMultiplier * smoothingMultiple
      next
    }

    # Use the provisional UMI threshold to estimate empty-droplet bounds.
    # These bounds anchor the debris-filtering step used for the refined
    # threshold estimate below.
    initialBounds <- findTrainingDataBoundsDefault(
      df_filtered,
      max_umis_empty,
      umiThresholdOverride = initialUmiThreshold,
      intronic_floor_fraction = intronic_floor_fraction,
      debris_pct_intronic_prior = debris_pct_intronic_prior,
      verbose = FALSE
    )

    if (!hasUsableEmptyBounds(initialBounds)) {
      resultDF <- addFailedSmoothingResult(
        resultDF = resultDF,
        smoothingMultiplier = smoothingMultiplier,
        umiThreshold = initialUmiThreshold
      )

      smoothingMultiplier <- smoothingMultiplier * smoothingMultiple
      next
    }

    # Second pass: remove likely high-UMI debris from the UMI threshold-learning
    # vector, but keep the full data set for final exemplar selection.
    thresholdLearningResult <- makeDebrisFilteredThresholdVector(
      df = df_filtered,
      bounds_empty = initialBounds$bounds_empty,
      debris_pct_intronic_prior = debris_pct_intronic_prior
    )

    umiThreshold <- determineUMIThreshold(
      thresholdLearningResult$x,
      method,
      denRange
    )

    if (is.na(umiThreshold)) {
      resultDF <- addFailedSmoothingResult(
        resultDF = resultDF,
        smoothingMultiplier = smoothingMultiplier,
        umiThreshold = NA_real_
      )

      smoothingMultiplier <- smoothingMultiplier * smoothingMultiple
      next
    }

    # Final pass: recompute empty and nucleus exemplar bounds using the refined
    # UMI threshold. The provisional bounds above are diagnostic only. Two
    # candidate debris floors are evaluated: the empty-derived default floor and
    # the biological debris prior.
    debrisFloorCandidates <- list(
      empty_derived = NULL,
      prior = debris_pct_intronic_prior
    )

    terminate_smoothing <- FALSE

    for (debrisFloorSource in names(debrisFloorCandidates)) {
      debrisFloorOverride <- debrisFloorCandidates[[debrisFloorSource]]

      bounds <- findTrainingDataBoundsDefault(
        df_filtered,
        max_umis_empty,
        umiThresholdOverride = umiThreshold,
        intronic_floor_fraction = intronic_floor_fraction,
        debris_pct_intronic_prior = debris_pct_intronic_prior,
        debris_intronic_floor_override = debrisFloorOverride,
        use2DTrainingRefinement = use2DTrainingRefinement,
        verbose = FALSE
      )

      if (!isUsableTrainingBoundsResult(bounds)) {
        resultDF <- addFailedSmoothingResult(
          resultDF = resultDF,
          smoothingMultiplier = smoothingMultiplier,
          umiThreshold = umiThreshold,
          debris_intronic_floor_source = debrisFloorSource
        )

        next
      }

      cell_features_labeled <- labelTrainingData(
        df_filtered,
        bounds$bounds_empty,
        bounds$bounds_non_empty,
        NULL,
        useCBRBInitialization = FALSE,
        training_empty_barcodes = bounds$training_empty_barcodes,
        training_nucleus_barcodes = bounds$training_nucleus_barcodes,
        verbose = FALSE
      )

      silhouetteResult <- calculate_silhouette(
        cell_features_labeled,
        showPlot = FALSE,
        verbose = FALSE
      )

      score <- silhouetteResult$mean_silhouette

      resultDF <- rbind(
        resultDF,
        makeSmoothingResultRow(
          smoothingMultiplier = smoothingMultiplier,
          umiThreshold = umiThreshold,
          debrisIntronicFloor = bounds$debris_intronic_floor,
          debrisIntronicFloorSource = debrisFloorSource,
          silhouetteScore = score
        )
      )

      if (!is.na(score)) {
        label_idx <- getTrainingLabelIndices(cell_features_labeled)

        candidate_results[[length(candidate_results) + 1]] <- list(
          resultIndex = nrow(resultDF),
          smoothingMultiplier = smoothingMultiplier,
          umiThreshold = umiThreshold,
          initial_umi_threshold = initialUmiThreshold,
          initial_empty_umi_upper_bound =
            thresholdLearningResult$initial_empty_umi_upper_bound,
          num_threshold_barcodes_removed =
            thresholdLearningResult$num_removed,
          debris_intronic_floor_source = debrisFloorSource,
          silhouette = score,
          bounds = bounds,
          debris_intronic_floor = bounds$debris_intronic_floor,
          empty_idx = label_idx$empty_idx,
          nucleus_idx = label_idx$nucleus_idx
        )

        if (score > best_silhouette) {
          best_silhouette <- score
        }

        # Early termination is disabled by default so all smoothing candidates
        # are evaluated. This allows later, more stable smoothing values to be
        # selected.
        if (!is.null(early_termination_threshold)) {
          threshold <- (1 - early_termination_threshold) * best_silhouette

          if (score < threshold) {
            log_info("Early termination: silhouette score dropped significantly")
            terminate_smoothing <- TRUE
            break
          }
        }
      }
    }

    if (terminate_smoothing) {
      break
    }
    # Advance smoothing after failed threshold/bounds attempts so a failed smoothing
    # value does not get retried repeatedly.
    smoothingMultiplier <- smoothingMultiplier * smoothingMultiple
  }

  candidate_results <- addCandidateStabilityScores(
    candidate_results,
    neighbor_window = stabilityNeighborWindow,
    max_stability_penalty = stabilityLambda
  )

  resultDF <- updateResultsWithCandidateScores(resultDF, candidate_results)
  best_candidate <- selectBestCandidateBySelectionScore(candidate_results)

  debris_result <- NULL
  if (!is.null(best_candidate)) {
    debris_result <- selectDebrisTrainingBounds(
      df = df_filtered,
      umiThreshold = best_candidate$umiThreshold,
      debris_intronic_floor = best_candidate$debris_intronic_floor,
      min_num = 10
    )
  }

  finalizeTrainingResults(
    df_filtered,
    best_candidate,
    resultDF,
    debris_result = debris_result
  )
}

#' Check whether final training bounds are usable for scoring
#'
#' @param bounds Result from `findTrainingDataBoundsDefault()`.
#'
#' @return Logical scalar.
#' @noRd
isUsableTrainingBoundsResult <- function(bounds) {
  if (is.null(bounds)) {
    return(FALSE)
  }

  if (is.null(bounds$debris_intronic_floor)) {
    return(FALSE)
  }

  if (is.null(bounds$bounds_empty) || is.null(bounds$bounds_non_empty)) {
    return(FALSE)
  }

  if (any(is.na(bounds$bounds_empty)) || any(is.na(bounds$bounds_non_empty))) {
    return(FALSE)
  }

  TRUE
}

#' Check whether initial empty bounds can be used for threshold refinement
#'
#' @param bounds Result from `findTrainingDataBoundsDefault()`.
#'
#' @return Logical scalar.
#' @noRd
hasUsableEmptyBounds <- function(bounds) {
  if (is.null(bounds) || is.null(bounds$bounds_empty)) {
    return(FALSE)
  }

  if (any(is.na(bounds$bounds_empty))) {
    return(FALSE)
  }

  TRUE
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

#' Make a debris-filtered UMI vector for threshold learning
#'
#' This helper uses provisional empty-droplet bounds to reduce the influence of
#' likely high-UMI debris during UMI threshold learning. A barcode is removed
#' from the threshold-learning vector only when it is both to the right of the
#' provisional empty-droplet UMI upper bound and below the debris pct_intronic
#' prior. This filtering is used only for threshold learning; final empty and
#' nucleus exemplar bounds are estimated from the full filtered data frame.
#'
#' @param df Data frame containing `num_transcripts` and `pct_intronic`.
#' @param bounds_empty Initial empty droplet bounds.
#' @param debris_pct_intronic_prior Numeric scalar. Prior expectation for the
#'     upper pct_intronic range of debris-like barcodes.
#'
#' @return A list containing the filtered UMI vector and diagnostic values.
#' @noRd
makeDebrisFilteredThresholdVector <- function(
  df,
  bounds_empty,
  debris_pct_intronic_prior = 0.25
) {
  x <- log10(df$num_transcripts + 1)

  drop <- x > bounds_empty$umi_upper_bound &
    df$pct_intronic < debris_pct_intronic_prior

  x_threshold <- x[!drop]

  if (length(x_threshold) == 0) {
    x_threshold <- x
    drop <- rep(FALSE, length(x))
  }

  list(
    x = x_threshold,
    initial_empty_umi_upper_bound = bounds_empty$umi_upper_bound,
    num_removed = sum(drop)
  )
}

getTrainingLabelIndices <- function(cell_features_labeled) {
  if (!"training_label_class" %in% colnames(cell_features_labeled)) {
    return(list(empty_idx = integer(0), nucleus_idx = integer(0)))
  }

  list(
    empty_idx = which(cell_features_labeled$training_label_class == "empty"),
    nucleus_idx = which(cell_features_labeled$training_label_class == "nucleus")
  )
}

#' Compute asymmetric containment between two index sets
#'
#' This function computes the fraction of entries in `x` that are also present
#' in `y`. It is asymmetric by design. If `y` is a strict superset of `x`, the
#' containment score is 1.
#'
#' @param x Integer vector.
#' @param y Integer vector.
#'
#' @return Numeric scalar containment score.
#' @noRd
containment <- function(x, y) {
  if (length(x) == 0) {
    return(NA_real_)
  }

  length(intersect(x, y)) / length(x)
}

#' Compute Jaccard similarity between two index sets
#'
#' This function computes the size of the intersection divided by the size of
#' the union. Unlike containment, Jaccard similarity penalizes expansions and
#' contractions of the selected exemplar set.
#'
#' @param x Integer vector.
#' @param y Integer vector.
#'
#' @return Numeric scalar Jaccard similarity.
#' @noRd
jaccard <- function(x, y) {
  union_length <- length(union(x, y))

  if (union_length == 0) {
    return(NA_real_)
  }

  length(intersect(x, y)) / union_length
}

#' Compute set similarity for candidate exemplar stability
#'
#' This function combines containment and Jaccard similarity by taking the
#' smaller value. Containment avoids over-penalizing small changes when the
#' current exemplar set is preserved nearby. Jaccard similarity penalizes large
#' expansions or contractions of the selected exemplar set.
#'
#' @param x Integer vector of exemplar indices for the current candidate.
#' @param y Integer vector of exemplar indices for a neighboring candidate.
#'
#' @return Numeric scalar similarity score.
#' @noRd
candidateSetSimilarity <- function(x, y) {
  containment_score <- containment(x, y)
  jaccard_score <- jaccard(x, y)

  min(containment_score, jaccard_score, na.rm = TRUE)
}


#' Get nearby candidates with the same debris floor source
#'
#' Candidate stability should be computed within a single debris-floor strategy.
#' This prevents candidates using the empty-derived floor from being compared
#' with candidates using the prior floor at the same or nearby smoothing values.
#'
#' @param candidate_results A list of candidate solution objects.
#' @param i Integer scalar. Index of the candidate of interest.
#' @param neighbor_window Integer scalar. Number of neighboring same-source
#'     candidates to include on each side when possible.
#' @param include_self Logical scalar. If `TRUE`, include candidate `i` in the
#'     returned window.
#'
#' @return Integer vector of candidate indices.
#' @noRd
getSourceMatchedCandidateWindow <- function(
  candidate_results,
  i,
  neighbor_window = 2,
  include_self = TRUE
) {
  source_i <- candidate_results[[i]]$debris_intronic_floor_source

  same_source_idx <- which(vapply(candidate_results, function(candidate) {
    identical(candidate$debris_intronic_floor_source, source_i)
  }, logical(1)))

  source_pos <- which(same_source_idx == i)

  if (length(source_pos) == 0) {
    return(integer(0))
  }

  local_pos <- getLocalWindowIndex(
    i = source_pos,
    n = length(same_source_idx),
    neighbor_window = neighbor_window
  )

  local_idx <- same_source_idx[local_pos]

  if (!include_self) {
    local_idx <- local_idx[local_idx != i]
  }

  local_idx
}


#' Compute local stability for a candidate solution
#'
#' This function estimates how stable one candidate exemplar selection is across
#' nearby smoothing values. The candidate is compared with neighboring
#' candidates within `neighbor_window` positions. Stability is computed
#' separately for empty droplet exemplars and nucleus exemplars, then summarized
#' by the smaller of the two values.
#'
#' The set similarity score combines containment and Jaccard similarity.
#' Containment asks whether the current exemplar set is preserved nearby, while
#' Jaccard similarity penalizes large expansions or contractions of the exemplar
#' set. This avoids treating a small candidate as fully stable when nearby
#' smoothing values select much broader regions.
#'
#' @param candidate_results A list of candidate solution objects. Each candidate
#'     is expected to contain `empty_idx` and `nucleus_idx` fields.
#' @param i Integer scalar. Index of the candidate whose stability should be
#'     computed.
#' @param neighbor_window Integer scalar. Number of neighboring candidates on
#'     each side to compare against.
#'
#' @return Numeric scalar stability score. Returns 1 when there are no
#'     neighboring candidates to compare. Returns `NA_real_` if stability cannot
#'     be computed.
#' @noRd
computeCandidateStability <- function(
  candidate_results,
  i,
  neighbor_window = 2
) {
  neighbor_idx <- getSourceMatchedCandidateWindow(
    candidate_results = candidate_results,
    i = i,
    neighbor_window = neighbor_window,
    include_self = FALSE
  )

  if (length(neighbor_idx) == 0) {
    return(1)
  }

  empty_i <- candidate_results[[i]]$empty_idx
  nucleus_i <- candidate_results[[i]]$nucleus_idx

  empty_scores <- vapply(neighbor_idx, function(j) {
    candidateSetSimilarity(empty_i, candidate_results[[j]]$empty_idx)
  }, numeric(1))

  nucleus_scores <- vapply(neighbor_idx, function(j) {
    candidateSetSimilarity(nucleus_i, candidate_results[[j]]$nucleus_idx)
  }, numeric(1))

  empty_stability <- mean(empty_scores, na.rm = TRUE)
  nucleus_stability <- mean(nucleus_scores, na.rm = TRUE)
  stability <- min(empty_stability, nucleus_stability, na.rm = TRUE)

  if (!is.finite(stability)) {
    return(NA_real_)
  }

  stability
}

#' Add stability-regularized selection scores to candidate solutions
#'
#' This function adds `stability`, `local_silhouette`, and `selection_score` to
#' each candidate solution. Stability measures how consistently the candidate's
#' selected exemplar sets are preserved across nearby smoothing values using
#' the same debris-floor source. The
#' local silhouette score is the median silhouette score among the candidate and
#' nearby smoothing values using the same debris-floor source. The final selection score subtracts a penalty for
#' instability from the local silhouette score.
#'
#' The penalty is bounded by `max_stability_penalty`. For example, if
#' `max_stability_penalty = 0.02`, then a completely unstable candidate can lose
#' at most 0.02 silhouette units. A fully stable candidate receives no penalty.
#'
#' @param candidate_results A list of candidate solution objects. Each candidate
#'     is expected to contain `silhouette`, `empty_idx`, and `nucleus_idx`.
#' @param neighbor_window Integer scalar. Number of nearby smoothing candidates
#'     on each side to use when estimating local stability and local silhouette.
#' @param max_stability_penalty Numeric scalar. Maximum amount subtracted from
#'     the local silhouette score when a candidate has zero stability.
#'
#' @return A list with the same structure as `candidate_results`. Each candidate
#'     has three additional fields: `stability`, `local_silhouette`, and
#'     `selection_score`.
#' @noRd
addCandidateStabilityScores <- function(
  candidate_results,
  neighbor_window = 2,
  max_stability_penalty = 0.02
) {
  if (length(candidate_results) == 0) {
    return(candidate_results)
  }

  for (i in seq_along(candidate_results)) {
    stability <- computeCandidateStability(
      candidate_results = candidate_results,
      i = i,
      neighbor_window = neighbor_window
    )

    local_silhouette <- computeLocalMedianSilhouette(
      candidate_results = candidate_results,
      i = i,
      neighbor_window = neighbor_window
    )

    candidate_results[[i]]$stability <- stability
    candidate_results[[i]]$local_silhouette <- local_silhouette
    candidate_results[[i]]$selection_score <-
      local_silhouette - max_stability_penalty * (1 - stability)
  }

  candidate_results
}


#' Add candidate stability scores to the per-iteration result table
#'
#' This function copies candidate-level stability metrics back into the
#' per-iteration result data frame. Rows without valid candidate solutions retain
#' missing values.
#'
#' @param results Data frame of per-iteration smoothing results.
#' @param candidate_results A list of candidate solution objects with
#'     `resultIndex`, `stability`, `local_silhouette`, and `selection_score`
#'     fields.
#'
#' @return The input `results` data frame with updated stability, local
#'     silhouette, and selection score columns.
#' @noRd
updateResultsWithCandidateScores <- function(results, candidate_results) {
  if (!"localSilhouette" %in% names(results)) {
    results$localSilhouette <- NA_real_
  }

  if (length(candidate_results) == 0) {
    return(results)
  }

  for (candidate in candidate_results) {
    i <- candidate$resultIndex
    results$stabilityScore[i] <- candidate$stability
    results$localSilhouette[i] <- candidate$local_silhouette
    results$selectionScore[i] <- candidate$selection_score
  }

  results
}

selectBestCandidateBySelectionScore <- function(candidate_results) {
  if (length(candidate_results) == 0) {
    return(NULL)
  }

  selection_scores <- vapply(candidate_results, function(candidate) {
    candidate$selection_score
  }, numeric(1))

  valid_idx <- which(!is.na(selection_scores))
  if (length(valid_idx) == 0) {
    return(NULL)
  }

  candidate_results[[valid_idx[which.max(selection_scores[valid_idx])]]]
}


finalizeTrainingResults <- function(
  df_filtered, best_candidate, results, debris_result = NULL
) {
  if (is.null(best_candidate)) {
    log_info("Unable to find good initialization. Returning Empty Bounds")
    return(list(
      best_silhouette = NA_real_,
      best_stability = NA_real_,
      best_selection_score = NA_real_,
      best_umi_threshold = NA_real_,
      initial_umi_threshold = NA_real_,
      initial_empty_umi_upper_bound = NA_real_,
      num_threshold_barcodes_removed = NA_real_,
      debris_intronic_floor = NA_real_,
      debris_intronic_floor_source = NA_character_,
      bounds_empty = NA,
      bounds_non_empty = NA,
      bounds_debris = makeEmptyTrainingBoundsDF(),
      training_empty_barcodes = character(0),
      training_nucleus_barcodes = character(0),
      training_debris_barcodes = character(0),
      resultDF = results,
      numEmpty = NA_real_,
      numNonEmpty = NA_real_,
      numDebris = 0
    ))
  }

  best_bounds <- best_candidate$bounds
  if (is.null(debris_result)) {
    debris_result <- makeEmptyDebrisTrainingBoundsResult()
  }

  cell_features_labeled <-
    labelTrainingData(df_filtered, best_bounds$bounds_empty,
      best_bounds$bounds_non_empty, NULL,
      useCBRBInitialization = FALSE,
      training_empty_barcodes = best_bounds$training_empty_barcodes,
      training_nucleus_barcodes = best_bounds$training_nucleus_barcodes,
      training_debris_barcodes = debris_result$training_debris_barcodes,
      bounds_debris = debris_result$bounds_debris,
      verbose = FALSE
    )

  numEmpty <- sum(cell_features_labeled$training_label_class == "empty",
    na.rm = TRUE
  )
  numNonEmpty <- sum(cell_features_labeled$training_label_class == "nucleus",
    na.rm = TRUE
  )
  numDebris <- sum(cell_features_labeled$training_label_class == "debris",
    na.rm = TRUE
  )

  return(list(
    best_silhouette = best_candidate$silhouette,
    best_stability = best_candidate$stability,
    best_selection_score = best_candidate$selection_score,
    best_umi_threshold = best_candidate$umiThreshold,
    initial_umi_threshold = best_candidate$initial_umi_threshold,
    initial_empty_umi_upper_bound =
      best_candidate$initial_empty_umi_upper_bound,
    num_threshold_barcodes_removed =
      best_candidate$num_threshold_barcodes_removed,
    debris_intronic_floor = best_candidate$debris_intronic_floor,
    debris_intronic_floor_source =
      best_candidate$debris_intronic_floor_source,
    bounds_empty = best_bounds$bounds_empty,
    bounds_non_empty = best_bounds$bounds_non_empty,
    bounds_debris = debris_result$bounds_debris,
    training_empty_barcodes = best_bounds$training_empty_barcodes,
    training_nucleus_barcodes = best_bounds$training_nucleus_barcodes,
    training_debris_barcodes = debris_result$training_debris_barcodes,
    resultDF = results,
    numEmpty = numEmpty,
    numNonEmpty = numNonEmpty,
    numDebris = numDebris
  ))
}


constrainDebrisBounds <- function(
  bounds,
  umiThreshold,
  debris_intronic_floor
) {
  bounds$umi_lower_bound <- max(bounds$umi_lower_bound, umiThreshold)
  bounds$intronic_lower_bound <- max(bounds$intronic_lower_bound, 0)
  bounds$intronic_upper_bound <- min(
    bounds$intronic_upper_bound,
    debris_intronic_floor
  )

  bounds
}

selectDebrisTrainingBounds <- function(
  df,
  umiThreshold,
  debris_intronic_floor,
  min_num = 10,
  pctDensity = 95,
  debris_start_quantile = 0.5,
  umi_upper_quantile = 0.995
) {
  empty_result <- makeEmptyDebrisTrainingBoundsResult()

  if (is.na(umiThreshold) || is.na(debris_intronic_floor)) {
    return(empty_result)
  }

  x <- log10(df$num_transcripts + 1)

  idx_debris_candidate <- which(
    x > umiThreshold &
      df$pct_intronic < debris_intronic_floor
  )

  if (length(idx_debris_candidate) < min_num) {
    return(empty_result)
  }

  debris_umi_floor <- as.numeric(stats::quantile(
    x[idx_debris_candidate],
    probs = debris_start_quantile,
    na.rm = TRUE
  ))

  idx_debris_for_density <- idx_debris_candidate[
    x[idx_debris_candidate] >= debris_umi_floor
  ]

  if (length(idx_debris_for_density) < min_num) {
    return(empty_result)
  }

  df_debris_density <- df[idx_debris_for_density, ]

  bounds_debris_intronic <- tryCatch(
    getHighestDensityIntervalsEnforcedSmoothing(
      df_debris_density,
      yAxisFeature = "pct_intronic",
      pctDensity = pctDensity,
      maxPeaksExpected = 1,
      showPlot = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(bounds_debris_intronic) || any(is.na(bounds_debris_intronic))) {
    return(empty_result)
  }

  bounds_debris_intronic <- constrainDebrisBounds(
    bounds = bounds_debris_intronic,
    umiThreshold = umiThreshold,
    debris_intronic_floor = debris_intronic_floor
  )

  df_debris_filtered <- df_debris_density[
    df_debris_density$pct_intronic >=
      bounds_debris_intronic$intronic_lower_bound &
      df_debris_density$pct_intronic <=
        bounds_debris_intronic$intronic_upper_bound,
  ]

  if (nrow(df_debris_filtered) < min_num) {
    return(empty_result)
  }

  bounds_debris_transcripts <- tryCatch(
    getHighestDensityIntervalsEnforcedSmoothing(
      df_debris_filtered,
      yAxisFeature = "pct_intronic",
      pctDensity = pctDensity,
      maxPeaksExpected = 1,
      showPlot = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(bounds_debris_transcripts) ||
    any(is.na(bounds_debris_transcripts))) {
    return(empty_result)
  }

  bounds_debris_transcripts <- constrainDebrisBounds(
    bounds = bounds_debris_transcripts,
    umiThreshold = umiThreshold,
    debris_intronic_floor = debris_intronic_floor
  )

  bounds_debris_transcripts$umi_upper_bound <- as.numeric(stats::quantile(
    log10(df_debris_filtered$num_transcripts + 1),
    umi_upper_quantile,
    na.rm = TRUE
  ))

  bounds_debris <- data.frame(
    umi_lower_bound = bounds_debris_transcripts$umi_lower_bound,
    umi_upper_bound = bounds_debris_transcripts$umi_upper_bound,
    intronic_lower_bound = bounds_debris_intronic$intronic_lower_bound,
    intronic_upper_bound = bounds_debris_intronic$intronic_upper_bound
  )

  if (bounds_debris$umi_lower_bound >= bounds_debris$umi_upper_bound ||
    bounds_debris$intronic_lower_bound >= bounds_debris$intronic_upper_bound) {
    return(empty_result)
  }

  selected <- log10(df_debris_density$num_transcripts + 1) >=
    bounds_debris$umi_lower_bound &
    log10(df_debris_density$num_transcripts + 1) <=
      bounds_debris$umi_upper_bound &
    df_debris_density$pct_intronic >= bounds_debris$intronic_lower_bound &
    df_debris_density$pct_intronic <= bounds_debris$intronic_upper_bound

  if (sum(selected) < min_num) {
    return(empty_result)
  }

  list(
    bounds_debris = bounds_debris,
    training_debris_barcodes = rownames(df_debris_density)[selected]
  )
}

makeEmptyDebrisTrainingBoundsResult <- function() {
  list(
    bounds_debris = makeEmptyTrainingBoundsDF(),
    training_debris_barcodes = character(0)
  )
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
#' This requires that there be at least 10 exemplar points in each partition
#' for the training bounds to be valid.  It's very hard to compute density
#' when the number of points is too low (especially 0 or 1) and that split
#' is unlikely to be a useful outcome.
#'
#' @param cell_features A dataframe containing the dataset with at least columns
#'   num_transcripts and pct_intronic.
#' @param max_umis_empty Exclude cell barcodes from analysis with fewer than
#'   this many UMIs. This is used to exclude very small barcodes that are noise
#'   and may confuse the initial turning point threshold.
#' @param umiThresholdOverride A numeric value to override the UMI threshold. If
#'   NULL, the function will find the UMI threshold using the turning point
#'   method. If not NULL, use an initial estimate provided by another means.
#' @param intronic_floor_fraction Fraction of the empty-droplet intronic lower
#'   bound used to filter the high-UMI partition before estimating nucleus
#'   bounds.
#' @param debris_pct_intronic_prior Numeric scalar. Prior expectation for the
#'   upper pct_intronic range of debris-like barcodes.
#' @param debris_intronic_floor_override Numeric scalar or `NULL`. If supplied,
#'   this value is used directly as the pct_intronic floor for candidate nuclei.
#'   If `NULL`, the floor is computed from the empty-derived intronic lower
#'   bound and the debris prior.
#' @param verbose Logical. if TRUE, prints verbose output to the log.
#'
#' @return A list with the training data bounds.
#' @noRd
findTrainingDataBoundsDefault <- function(
  cell_features, max_umis_empty = 50,
  umiThresholdOverride = NULL,
  intronic_floor_fraction = 1,
  debris_pct_intronic_prior = 0.25,
  debris_intronic_floor_override = NULL,
  use2DTrainingRefinement = FALSE,
  verbose = TRUE
) {
  df <- cell_features[cell_features$num_transcripts >= max_umis_empty, ]

  umiThreshold <- if (!is.null(umiThresholdOverride)) {
    umiThresholdOverride
  } else {
    PitAfterHighestPeakWithGridSearch(log10(df$num_transcripts + 1))
  }

  if (verbose) {
    log_info("Initial UMI Threshold [", round(umiThreshold, 3), "]")
  }

  df_empty <- df[log10(df$num_transcripts + 1) < umiThreshold, ]
  bounds_empty <- getEmptyCellsByDensity(df_empty, "pct_intronic",
    pctDensity = 75, showPlot = FALSE, verbose = verbose
  )

  early_exit <- checkEarlyExit(df_empty, bounds_empty, umiThreshold,
    min_num = 10, verbose
  )
  if (!is.null(early_exit)) {
    return(early_exit)
  }

  # Keep the empty exemplar bounds consistent with the UMI partition.
  # Density smoothing can extend the interval slightly past the threshold.
  bounds_empty$umi_upper_bound <- min(
    bounds_empty$umi_upper_bound,
    umiThreshold
  )

  # Capture the pct_intronic floor used for candidate nuclei.
  non_empty_filter_result <- filterNonEmptyPartitionByEmptyIntronicBound(
    df = df,
    umiThreshold = umiThreshold,
    bounds_empty = bounds_empty,
    intronic_floor_fraction = intronic_floor_fraction,
    debris_pct_intronic_prior = debris_pct_intronic_prior,
    debris_intronic_floor_override = debris_intronic_floor_override,
    verbose = verbose
  )

  df_non_empty <- non_empty_filter_result$df_non_empty
  debris_intronic_floor <- non_empty_filter_result$debris_intronic_floor

  early_exit <- checkEarlyExit(
    df_non_empty, bounds_empty,
    umiThreshold,
    min_num = 10, verbose
  )
  if (!is.null(early_exit)) {
    early_exit$debris_intronic_floor <- debris_intronic_floor
    return(early_exit)
  }


  bounds_non_empty_intronic <- getHighestDensityIntervalsEnforcedSmoothing(
    df_non_empty,
    yAxisFeature = "pct_intronic", pctDensity = 75,
    maxPeaksExpected = 1, showPlot = FALSE
  )
  df_filtered <- df_non_empty[
    df_non_empty$pct_intronic >=
      bounds_non_empty_intronic$intronic_lower_bound &
      df_non_empty$pct_intronic <=
        bounds_non_empty_intronic$intronic_upper_bound,
  ]

  bounds_non_empty_transcripts <- getHighestDensityIntervalsEnforcedSmoothing(
    df_filtered,
    yAxisFeature = "pct_intronic",
    pctDensity = 75, showPlot = FALSE
  )

  # Extend the upper UMI bound to include the high-UMI tail within the selected
  # nucleus intronic interval. High-UMI barcodes in this interval are likely
  # nuclei even if they are too sparse to belong to the highest-density UMI
  # interval.
  bounds_non_empty_transcripts$umi_upper_bound <-
    as.numeric(stats::quantile(log10(df_filtered$num_transcripts + 1), 0.995))

  bounds_non_empty <- data.frame(
    umi_lower_bound = bounds_non_empty_transcripts$umi_lower_bound,
    umi_upper_bound = bounds_non_empty_transcripts$umi_upper_bound,
    intronic_lower_bound = bounds_non_empty_intronic$intronic_lower_bound,
    intronic_upper_bound = bounds_non_empty_intronic$intronic_upper_bound
  )

  # Keep the nucleus exemplar bounds consistent with the UMI partition.
  # Density smoothing can extend the interval slightly below the threshold.
  bounds_non_empty$umi_lower_bound <- max(
    bounds_non_empty$umi_lower_bound,
    umiThreshold
  )

  training_empty_barcodes <- NULL
  training_nucleus_barcodes <- NULL
  bounds_empty_2d <- NULL
  bounds_non_empty_2d <- NULL

  if (use2DTrainingRefinement) {
    empty_refined <- refineTrainingBoundsWith2DComponent(
      df = df_empty,
      bounds = bounds_empty,
      pctDensity = 75
    )

    nucleus_refined <- refineTrainingBoundsWith2DComponent(
      df = df_non_empty,
      bounds = bounds_non_empty,
      pctDensity = 75
    )

    if (!isUsable2DRefinement(empty_refined) ||
      !isUsable2DRefinement(nucleus_refined)) {
      return(makeEmptyTrainingBoundsResult(
        bounds_empty = bounds_empty,
        debris_intronic_floor = debris_intronic_floor
      ))
    }

    training_empty_barcodes <- empty_refined$selected_barcodes
    training_nucleus_barcodes <- extendNucleus2DSelectionToHighUMITail(
      df = df_non_empty,
      selected_barcodes = nucleus_refined$selected_barcodes
    )

    bounds_empty_2d <- empty_refined$bounds
    bounds_non_empty_2d <- getBoundsFromSelectedBarcodes(
      df = df_non_empty,
      selected_barcodes = training_nucleus_barcodes
    )
  }

  return(list(
    bounds_empty = bounds_empty,
    bounds_non_empty = bounds_non_empty,
    bounds_empty_2d = bounds_empty_2d,
    bounds_non_empty_2d = bounds_non_empty_2d,
    training_empty_barcodes = training_empty_barcodes,
    training_nucleus_barcodes = training_nucleus_barcodes,
    debris_intronic_floor = debris_intronic_floor
  ))
}

extendNucleus2DSelectionToHighUMITail <- function(
  df,
  selected_barcodes,
  high_umi_quantile = 0.75,
  intronic_quantiles = c(0.05, 0.95),
  umi_upper_quantile = 0.995
) {
  if (length(selected_barcodes) == 0) {
    return(character(0))
  }

  selected <- rownames(df) %in% selected_barcodes

  if (sum(selected) == 0) {
    return(character(0))
  }

  x <- log10(df$num_transcripts + 1)
  y <- df$pct_intronic

  selected_x <- x[selected]
  selected_y <- y[selected]

  high_umi_cutoff <- as.numeric(stats::quantile(
    selected_x,
    high_umi_quantile,
    na.rm = TRUE
  ))

  high_umi_selected <- selected & x >= high_umi_cutoff

  if (sum(high_umi_selected) < 10) {
    high_umi_selected <- selected
  }

  intronic_bounds <- as.numeric(stats::quantile(
    y[high_umi_selected],
    intronic_quantiles,
    na.rm = TRUE
  ))

  intronic_lower <- intronic_bounds[1]
  intronic_upper <- intronic_bounds[2]

  umi_lower <- min(selected_x, na.rm = TRUE)

  in_intronic_band <- y >= intronic_lower & y <= intronic_upper

  if (sum(in_intronic_band) == 0) {
    return(selected_barcodes)
  }

  umi_upper <- as.numeric(stats::quantile(
    x[in_intronic_band],
    umi_upper_quantile,
    na.rm = TRUE
  ))

  extended <- x >= umi_lower &
    x <= umi_upper &
    in_intronic_band

  union(selected_barcodes, rownames(df)[extended])
}

getBoundsFromSelectedBarcodes <- function(df, selected_barcodes) {
  selected <- rownames(df) %in% selected_barcodes

  if (sum(selected) == 0) {
    return(makeEmptyTrainingBoundsDF())
  }

  x <- log10(df$num_transcripts + 1)
  y <- df$pct_intronic

  data.frame(
    umi_lower_bound = min(x[selected], na.rm = TRUE),
    umi_upper_bound = max(x[selected], na.rm = TRUE),
    intronic_lower_bound = min(y[selected], na.rm = TRUE),
    intronic_upper_bound = max(y[selected], na.rm = TRUE)
  )
}


#' Refine rectangular exemplar bounds with a 2D HDR component
#'
#' This function estimates a 2D highest-density region in log10 UMI and
#' pct_intronic space. If the HDR contains multiple connected components, the
#' component with the largest overlap with the supplied rectangular bounds is
#' selected. The selected cell barcodes are returned using row names.
#'
#' @param df Data frame containing `num_transcripts` and `pct_intronic`.
#' @param bounds Rectangular bounds used as the anchor for choosing the 2D
#'     connected component.
#' @param pctDensity Numeric scalar. Probability mass used for the 2D HDR.
#' @param kde_package Character scalar passed to `hdrcde::hdr.2d()`.
#'
#' @return A list containing refined bounds, selected barcodes, and diagnostics.
#' @noRd
refineTrainingBoundsWith2DComponent <- function(
  df,
  bounds,
  pctDensity = 75,
  kde_package = "ash"
) {
  empty_result <- makeEmpty2DRefinementResult()

  if (is.null(rownames(df)) || any(rownames(df) == "")) {
    return(empty_result)
  }

  x <- log10(df$num_transcripts + 1)
  y <- df$pct_intronic

  keep <- is.finite(x) & is.finite(y)
  x_in <- x[keep]
  y_in <- y[keep]

  if (length(x_in) < 10) {
    return(empty_result)
  }

  hdr <- tryCatch(
    hdrcde::hdr.2d(
      x = x_in,
      y = y_in,
      prob = pctDensity / 100,
      kde.package = kde_package
    ),
    error = function(e) NULL
  )

  if (is.null(hdr)) {
    return(empty_result)
  }

  density_threshold <- as.numeric(hdr$falpha[1])
  density_mask <- hdr$den$z >= density_threshold
  component_grid <- labelConnectedGridComponents(density_mask)

  point_component <- assignPointsToGridComponents(
    x = x_in,
    y = y_in,
    grid_x = hdr$den$x,
    grid_y = hdr$den$y,
    component_grid = component_grid
  )

  in_anchor <- x_in >= bounds$umi_lower_bound &
    x_in <= bounds$umi_upper_bound &
    y_in >= bounds$intronic_lower_bound &
    y_in <= bounds$intronic_upper_bound

  component_ids <- sort(unique(point_component))
  component_ids <- component_ids[!is.na(component_ids) & component_ids > 0]

  if (length(component_ids) == 0 || sum(in_anchor) == 0) {
    empty_result$hdr <- hdr
    return(empty_result)
  }

  component_summary <- summarize2DComponents(
    component_ids = component_ids,
    point_component = point_component,
    in_anchor = in_anchor,
    x_in = x_in,
    y_in = y_in
  )

  if (nrow(component_summary) == 0 ||
    max(component_summary$n_overlap, na.rm = TRUE) == 0) {
    empty_result$hdr <- hdr
    empty_result$component_summary <- component_summary
    return(empty_result)
  }

  best_component <- component_summary$component_id[
    which.max(component_summary$n_overlap)
  ]

  selected_in <- point_component == best_component

  if (!is.null(hdr$fxy)) {
    selected_in <- selected_in & hdr$fxy >= density_threshold
  }

  selected <- rep(FALSE, nrow(df))
  selected[keep] <- selected_in

  selected_barcodes <- rownames(df)[selected]

  if (length(selected_barcodes) == 0) {
    empty_result$hdr <- hdr
    empty_result$component_summary <- component_summary
    return(empty_result)
  }

  list(
    bounds = data.frame(
      umi_lower_bound = min(x[selected]),
      umi_upper_bound = max(x[selected]),
      intronic_lower_bound = min(y[selected]),
      intronic_upper_bound = max(y[selected])
    ),
    selected = selected,
    selected_barcodes = selected_barcodes,
    hdr = hdr,
    component_summary = component_summary,
    selected_component = best_component,
    density_threshold = density_threshold,
    pctDensity = pctDensity,
    kde_package = kde_package,
    anchor_bounds = bounds
  )
}


summarize2DComponents <- function(
  component_ids,
  point_component,
  in_anchor,
  x_in,
  y_in
) {
  do.call(
    rbind,
    lapply(component_ids, function(component_id) {
      in_component <- point_component == component_id
      n_component <- sum(in_component)
      n_overlap <- sum(in_component & in_anchor)

      data.frame(
        component_id = component_id,
        n_component = n_component,
        component_fraction = n_component / length(x_in),
        n_overlap = n_overlap,
        overlap_fraction = n_overlap / n_component,
        anchor_recall = n_overlap / sum(in_anchor),
        median_umi = stats::median(x_in[in_component]),
        median_intronic = stats::median(y_in[in_component])
      )
    })
  )
}

labelConnectedGridComponents <- function(mask) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  labels <- matrix(0L, nrow = nr, ncol = nc)
  current_label <- 0L

  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      if (!mask[i, j] || labels[i, j] != 0L) {
        next
      }

      current_label <- current_label + 1L
      queue_i <- i
      queue_j <- j
      labels[i, j] <- current_label

      while (length(queue_i) > 0) {
        qi <- queue_i[1]
        qj <- queue_j[1]
        queue_i <- queue_i[-1]
        queue_j <- queue_j[-1]

        neighbors <- rbind(
          c(qi - 1L, qj),
          c(qi + 1L, qj),
          c(qi, qj - 1L),
          c(qi, qj + 1L)
        )

        for (k in seq_len(nrow(neighbors))) {
          ni <- neighbors[k, 1]
          nj <- neighbors[k, 2]

          if (ni < 1L || ni > nr || nj < 1L || nj > nc) {
            next
          }

          if (!mask[ni, nj] || labels[ni, nj] != 0L) {
            next
          }

          labels[ni, nj] <- current_label
          queue_i <- c(queue_i, ni)
          queue_j <- c(queue_j, nj)
        }
      }
    }
  }

  labels
}

assignPointsToGridComponents <- function(
  x,
  y,
  grid_x,
  grid_y,
  component_grid
) {
  x_idx <- findInterval(x, grid_x, all.inside = TRUE)
  y_idx <- findInterval(y, grid_y, all.inside = TRUE)

  component_grid[cbind(x_idx, y_idx)]
}

makeEmpty2DRefinementResult <- function() {
  list(
    bounds = makeEmptyTrainingBoundsDF(),
    selected = logical(0),
    selected_barcodes = character(0),
    hdr = NULL,
    component_summary = NULL,
    selected_component = NA_integer_,
    density_threshold = NA_real_,
    pctDensity = NA_real_,
    kde_package = NA_character_,
    anchor_bounds = NULL
  )
}

isUsable2DRefinement <- function(result) {
  if (is.null(result) || length(result$selected_barcodes) == 0) {
    return(FALSE)
  }

  if (is.null(result$bounds) || any(is.na(result$bounds))) {
    return(FALSE)
  }

  TRUE
}

makeEmptyTrainingBoundsDF <- function() {
  data.frame(
    umi_lower_bound = NA_real_,
    umi_upper_bound = NA_real_,
    intronic_lower_bound = NA_real_,
    intronic_upper_bound = NA_real_
  )
}

makeEmptyTrainingBoundsResult <- function(
  bounds_empty = makeEmptyTrainingBoundsDF(),
  debris_intronic_floor = NA_real_
) {
  list(
    bounds_empty = bounds_empty,
    bounds_non_empty = makeEmptyTrainingBoundsDF(),
    bounds_empty_2d = makeEmptyTrainingBoundsDF(),
    bounds_non_empty_2d = makeEmptyTrainingBoundsDF(),
    bounds_debris = makeEmptyTrainingBoundsDF(),
    training_empty_barcodes = character(0),
    training_nucleus_barcodes = character(0),
    training_debris_barcodes = character(0),
    debris_intronic_floor = debris_intronic_floor
  )
}


filterNonEmptyPartitionByEmptyIntronicBound <- function(
  df,
  umiThreshold,
  bounds_empty,
  intronic_floor_fraction = 1,
  debris_pct_intronic_prior = 0.25,
  debris_intronic_floor_override = NULL,
  verbose = FALSE
) {
  df_high_umi <- df[log10(df$num_transcripts + 1) > umiThreshold, ]

  observed_intronic_floor <- intronic_floor_fraction *
    bounds_empty$intronic_lower_bound

  debris_intronic_floor <- if (!is.null(debris_intronic_floor_override)) {
    debris_intronic_floor_override
  } else {
    min(
      observed_intronic_floor,
      debris_pct_intronic_prior
    )
  }

  if (verbose) {
    logger::log_info(
      paste0(
        "Debris pct intronic prior [",
        round(debris_pct_intronic_prior, 3),
        "], observed empty-derived intronic floor [",
        round(observed_intronic_floor, 3),
        "], applied debris intronic floor [",
        round(debris_intronic_floor, 3),
        "]"
      )
    )
  }

  list(
    df_non_empty = df_high_umi[
      df_high_umi$pct_intronic >= debris_intronic_floor,
    ],
    debris_intronic_floor = debris_intronic_floor,
    observed_intronic_floor = observed_intronic_floor
  )
}

checkEarlyExit <- function(df, bounds_empty, umiThreshold, min_num = 10, verbose) {
  if (nrow(df) < min_num || is.null(bounds_empty)) {
    if (verbose) log_warn("Minimal data left after selecting barcodes for one partition.")
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

# If there are multiple peaks, only keep the largest peak.
getEmptyCellsByDensity <- function(
  cell_features, yAxisFeature = "pct_intronic",
  pctDensity = 75, showPlot = FALSE, verbose = FALSE
) {
  probList <- unique(c(25, 50, 75, 90, 95, 99, pctDensity))
  probList <- probList[probList <= pctDensity]

  # Do this by density.  This distribution is bimodal (empty/nuclei.)
  x <- log10(cell_features$num_transcripts + 1)
  zHDR <- getBoundsByDensity(x,
    probList = probList, pctDensity = pctDensity,
    maxPeaksExpected = 2, showPlot = showPlot
  )

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
      msg <- paste(
        "Unable to select the empty cell barcode highest",
        "density interval."
      )
      log_warn(msg)
    }
    return(NULL)
  }
  umiBounds <- zHDR_intervals[max_idx, ]

  # select data from this interval and find the %intronic.
  df <- cell_features[
    log10(cell_features$num_transcripts + 1) >= umiBounds[1] &
      log10(cell_features$num_transcripts + 1) <= umiBounds[2],
  ]
  pctDensity <- 95
  probList <- unique(c(25, 50, 75, 90, 95, 99, pctDensity))
  probList <- probList[probList <= pctDensity]
  yBounds <- getBoundsByDensity(df[[yAxisFeature]], probList, pctDensity,
    maxPeaksExpected = 1, showPlot = showPlot
  )$intervals
  result <- data.frame(
    umi_lower_bound = umiBounds[1],
    umi_upper_bound = umiBounds[2], intronic_lower_bound = yBounds[1],
    intronic_upper_bound = yBounds[2]
  )
  return(result)
}

#' Compute local median silhouette for a candidate solution
#'
#' This function summarizes the silhouette score for a candidate and its nearby
#' smoothing neighbors. Using the local median silhouette reduces the influence
#' of isolated silhouette spikes that are not supported by neighboring smoothing
#' values.
#'
#' @param candidate_results A list of candidate solution objects. Each candidate
#'     is expected to contain a `silhouette` field.
#' @param i Integer scalar. Index of the candidate whose local silhouette should
#'     be computed.
#' @param neighbor_window Integer scalar. Number of neighboring candidates on
#'     each side to include.
#'
#' @return Numeric scalar local median silhouette score.
#' @noRd
computeLocalMedianSilhouette <- function(
  candidate_results,
  i,
  neighbor_window = 2
) {
  local_idx <- getSourceMatchedCandidateWindow(
    candidate_results = candidate_results,
    i = i,
    neighbor_window = neighbor_window,
    include_self = TRUE
  )

  silhouette_values <- vapply(local_idx, function(j) {
    candidate_results[[j]]$silhouette
  }, numeric(1))

  stats::median(silhouette_values, na.rm = TRUE)
}

#' Create an empty smoothing-result data frame
#'
#' @return An empty data frame with the columns used for iterative smoothing
#'     diagnostics.
#' @noRd
makeEmptySmoothingResultDF <- function() {
  data.frame(
    smoothingMultiplier = numeric(0),
    umiThreshold = numeric(0),
    debrisIntronicFloor = numeric(0),
    debrisIntronicFloorSource = character(0),
    silhouetteScore = numeric(0),
    localSilhouette = numeric(0),
    stabilityScore = numeric(0),
    selectionScore = numeric(0)
  )
}

#' Create one smoothing-result row
#'
#' @param smoothingMultiplier Numeric scalar smoothing multiplier.
#' @param umiThreshold Numeric scalar UMI threshold.
#' @param silhouetteScore Numeric scalar silhouette score.
#' @param debrisIntronicFloor Numeric scalar debris/nucleus intronic floor.
#' @param debrisIntronicFloorSource Character scalar describing the floor source.
#' @param localSilhouette Numeric scalar local silhouette score.
#' @param stabilityScore Numeric scalar stability score.
#' @param selectionScore Numeric scalar stability-regularized selection score.
#'
#' @return A one-row data frame with fixed diagnostic column names.
#' @noRd
makeSmoothingResultRow <- function(
  smoothingMultiplier,
  umiThreshold,
  silhouetteScore,
  debrisIntronicFloor = NA_real_,
  debrisIntronicFloorSource = NA_character_,
  localSilhouette = NA_real_,
  stabilityScore = NA_real_,
  selectionScore = NA_real_
) {
  debrisIntronicFloor <- scalarOrDefault(debrisIntronicFloor, NA_real_)
  debrisIntronicFloorSource <- scalarOrDefault(
    debrisIntronicFloorSource,
    NA_character_
  )

  data.frame(
    smoothingMultiplier = smoothingMultiplier,
    umiThreshold = umiThreshold,
    debrisIntronicFloor = debrisIntronicFloor,
    debrisIntronicFloorSource = debrisIntronicFloorSource,
    silhouetteScore = silhouetteScore,
    localSilhouette = localSilhouette,
    stabilityScore = stabilityScore,
    selectionScore = selectionScore
  )
}

#' Return a scalar value or a default when the input is empty
#'
#' @param x Input value.
#' @param default Default value used when `x` is `NULL` or length zero.
#'
#' @return A scalar value.
#' @noRd
scalarOrDefault <- function(x, default) {
  if (is.null(x) || length(x) == 0) {
    return(default)
  }

  x[1]
}

#' Add a failed smoothing-result row
#'
#' @param resultDF Data frame of per-iteration smoothing results.
#' @param smoothingMultiplier Numeric scalar smoothing multiplier.
#' @param umiThreshold Numeric scalar UMI threshold.
#' @param debris_intronic_floor_source Character scalar describing the floor source.
#'
#' @return `resultDF` with one additional row containing missing scores.
#' @noRd
addFailedSmoothingResult <- function(
  resultDF,
  smoothingMultiplier,
  umiThreshold,
  debris_intronic_floor_source = NA_character_
) {
  rbind(
    resultDF,
    makeSmoothingResultRow(
      smoothingMultiplier = smoothingMultiplier,
      umiThreshold = umiThreshold,
      debrisIntronicFloorSource = debris_intronic_floor_source,
      silhouetteScore = NA_real_
    )
  )
}
