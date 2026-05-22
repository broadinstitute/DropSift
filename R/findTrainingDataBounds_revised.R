####################
# HIGH LEVEL BOUNDS FUNCTIONS
#####################

#' Find exemplars for nuclei and empty droplets
#'
#' This function selects exemplar barcodes for empty droplets and nuclei. These
#' exemplars are used downstream to train the SVM classifier. When CellBender
#' features are not used, the algorithm first estimates a UMI threshold that
#' separates the likely-empty and likely-nucleus barcode partitions. It then
#' selects high-density exemplar regions within each partition using log10 UMI
#' counts and pct_intronic.
#'
#' In the non-CellBender path, low-intronic debris can sometimes fall in the
#' high-UMI partition and distort nucleus exemplar selection. To reduce this
#' failure mode, the algorithm applies an intronic floor before estimating the
#' nucleus exemplar bounds. The observed floor is computed from the empty droplet
#' exemplar bounds:
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
#' @param cell_features The cell features data frame.
#' @param max_umis_empty Exclude cell barcodes from analysis with fewer than
#'   this many UMIs. This threshold should exclude noisy barcodes, but not
#'   exclude the empty droplet cloud.
#' @param useCellBenderFeatures When true, the CellBender remove-background
#'   feature frac_contamination is used for cell selection.
#' @param forceTwoClusterSolution When true, the function will attempt to find a
#'   solution with two clusters. This may be useful when the data is overloaded
#'   and breaks the normal assumptions, but may find suboptimal solutions for
#'   other data sets.
#' @param intronic_floor_fraction Numeric scalar. Multiplier applied to the
#'   empty droplet intronic lower bound when computing the observed intronic
#'   floor used to remove low-intronic debris from the candidate nucleus
#'   partition.
#' @param debris_pct_intronic_prior Numeric scalar. Prior expectation for the
#'   upper pct_intronic range of debris-like barcodes. When filtering the
#'   high-UMI candidate nucleus partition, the algorithm uses the smaller of
#'   this value and the empty-derived intronic lower bound learned from the
#'   empty droplet density as the intronic floor. This prevents a high-intronic
#'   empty droplet cluster from imposing an unrealistically high debris-removal
#'   threshold.  Operationally, candidate nuclei must have pct_intronic
#'   greater than or equal to the smaller of this prior and the adjusted
#'   empty-derived intronic lower bound.
#' @param verbose Print verbose output to log.
#'
#' @return A list containing the bounds for the empty droplet and nuclei
#'   exemplars.
#' @import logger
#' @noRd
findTrainingDataBounds <- function(
  cell_features, max_umis_empty = 50,
  useCellBenderFeatures = TRUE, forceTwoClusterSolution = FALSE,
  intronic_floor_fraction = 0.9, debris_pct_intronic_prior = 0.25,
  verbose = FALSE
) {
  # If using CellBender features, use the specific CBRB method.
  if (useCellBenderFeatures) {
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
      debris_pct_intronic_prior, verbose
    ))
  } else {
    # Otherwise, use a more general approach with empty droplet fraction
    # constraints.
    return(findDefaultSolution(
      cell_features, max_umis_empty,
      max_umis_empty_off, intronic_floor_fraction,
      debris_pct_intronic_prior, verbose
    ))
  }
}

findTwoClusterSolution <- function(
  cell_features, max_umis_empty,
  max_umis_empty_off, intronic_floor_fraction,
  debris_pct_intronic_prior, verbose
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
  debris_pct_intronic_prior, verbose
) {
  # The default approach finds the pit after the highest peak.
  logSelectionProcess("Pit After Highest Peak", max_umis_empty, verbose)
  default <- findTrainingDataBoundsDefaultIterative(
    cell_features,
    max_umis_empty = max_umis_empty,
    intronic_floor_fraction = intronic_floor_fraction,
    debris_pct_intronic_prior = debris_pct_intronic_prior,
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

#' Find the training data bounds using an iterative approach
#'
#' This function scans smoothing values, estimates a UMI threshold for each
#' smoothing value, selects candidate exemplar bounds, and scores each candidate
#' using the silhouette score. Candidate solutions are then assigned local
#' stability scores based on whether their selected exemplars are preserved
#' across nearby smoothing values. The final candidate is selected using a
#' stability-regularized score.
#'
#' @param cell_features A data frame containing barcode-level features.
#' @param max_umis_empty Exclude cell barcodes from analysis with fewer than
#'     this many UMIs before estimating the initial UMI threshold.
#' @param method Method used to select the UMI threshold.
#' @param smoothingMultiple Multiplicative factor used to increase smoothing
#'     across iterations.
#' @param max_iterations Maximum number of smoothing iterations.
#' @param early_termination_threshold Fractional silhouette drop used for early
#'     stopping.
#' @param stabilityLambda Maximum silhouette penalty applied to a completely
#'     unstable candidate.
#' @param stabilityNeighborWindow Number of nearby candidates on each side used
#'     to estimate local stability.
#' @param intronic_floor_fraction Fraction of the empty-droplet intronic lower
#'     bound used to filter the high-UMI partition before estimating nucleus
#'     bounds.
#' @param debris_pct_intronic_prior Numeric scalar. Prior expectation for the
#'     upper pct_intronic range of debris-like barcodes.
#' @param verbose Logical scalar. If `TRUE`, log progress.
#'
#' @return A list containing selected bounds, selection metrics, and the
#'     per-iteration diagnostic data frame.
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
    verbose = TRUE) {
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

    if (is.null(initialBounds)) {
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
    # UMI threshold. The provisional bounds above are diagnostic only.
    bounds <- findTrainingDataBoundsDefault(
      df_filtered,
      max_umis_empty,
      umiThresholdOverride = umiThreshold,
      intronic_floor_fraction = intronic_floor_fraction,
      debris_pct_intronic_prior = debris_pct_intronic_prior,
      verbose = FALSE
    )

    if (is.null(bounds)) {
      resultDF <- addFailedSmoothingResult(
        resultDF = resultDF,
        smoothingMultiplier = smoothingMultiplier,
        umiThreshold = umiThreshold
      )

      smoothingMultiplier <- smoothingMultiplier * smoothingMultiple
      next
    }

    cell_features_labeled <- labelTrainingData(
      df_filtered,
      bounds$bounds_empty,
      bounds$bounds_non_empty,
      NULL,
      useCellBenderFeatures = FALSE,
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
        silhouette = score,
        bounds = bounds,
        debris_intronic_floor = bounds$debris_intronic_floor,
        empty_idx = label_idx$empty_idx,
        nucleus_idx = label_idx$nucleus_idx
      )

      if (score > best_silhouette) {
        best_silhouette <- score
      }

      threshold <- (1 - early_termination_threshold) * best_silhouette

      # Early termination is disabled by default so all smoothing candidates are
      # evaluated. This makes candidate selection less dependent on the first local
      # silhouette drop and allows later, more stable smoothing values to be selected.
      if (!is.null(early_termination_threshold)) {
        threshold <- (1 - early_termination_threshold) * best_silhouette

        if (score < threshold) {
          log_info("Early termination: silhouette score dropped significantly")
          break
        }
      }
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

  finalizeTrainingResults(
    df_filtered,
    best_candidate,
    resultDF
  )
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
    debris_pct_intronic_prior = 0.25) {

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
  is_labeled <- !is.na(cell_features_labeled$training_label_is_cell)

  list(
    empty_idx = which(is_labeled &
      !cell_features_labeled$training_label_is_cell),
    nucleus_idx = which(is_labeled &
      cell_features_labeled$training_label_is_cell)
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
  if (length(candidate_results) < 2) {
    return(1)
  }

  neighbor_idx <- getLocalWindowIndex(
    i = i,
    n = length(candidate_results),
    neighbor_window = neighbor_window
  )

  neighbor_idx <- neighbor_idx[neighbor_idx != i]

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
#' selected exemplar sets are preserved across nearby smoothing values. The
#' local silhouette score is the median silhouette score among the candidate and
#' nearby smoothing values. The final selection score subtracts a penalty for
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
  df_filtered, best_candidate, results
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
      bounds_empty = NA,
      bounds_non_empty = NA,
      resultDF = results,
      numEmpty = NA_real_,
      numNonEmpty = NA_real_
    ))
  }

  best_bounds <- best_candidate$bounds
  cell_features_labeled <-
    labelTrainingData(df_filtered, best_bounds$bounds_empty,
      best_bounds$bounds_non_empty, NULL,
      useCellBenderFeatures = FALSE, verbose = FALSE
    )

  numEmpty <- sum(!is.na(cell_features_labeled$training_label_is_cell) &
    !cell_features_labeled$training_label_is_cell)
  numNonEmpty <- sum(!is.na(cell_features_labeled$training_label_is_cell) &
    cell_features_labeled$training_label_is_cell)

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
    bounds_empty = best_bounds$bounds_empty,
    bounds_non_empty = best_bounds$bounds_non_empty,
    resultDF = results,
    numEmpty = numEmpty,
    numNonEmpty = numNonEmpty
  ))
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
#' @param debris_pct_intronic_prior Numeric scalar. Maximum empty-droplet
#'   intronic lower bound that is interpreted as a low-intronic debris boundary.
#'   If the computed floor is above this value, the floor is not applied.
#' @param verbose Logical. if TRUE, prints verbose output to the log.
#'
#' @return A list with the training data bounds.
#' @noRd
findTrainingDataBoundsDefault <- function(
  cell_features, max_umis_empty = 50,
  umiThresholdOverride = NULL,
  intronic_floor_fraction = 1,
  debris_pct_intronic_prior = 0.25,
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

  # TODO: added a minimum number of cell barcodes that must be present
  # in each partition when checking if this result is valid,
  # instead of checking for 0 cell barcodes in the partition.
  early_exit <- checkEarlyExit(df_empty, bounds_empty, umiThreshold,
    min_num = 10, verbose
  )
  if (!is.null(early_exit)) {
    return(early_exit)
  }

  # Capture the pct_intronic floor used for candidate nuclei.
  non_empty_filter_result <- filterNonEmptyPartitionByEmptyIntronicBound(
    df = df,
    umiThreshold = umiThreshold,
    bounds_empty = bounds_empty,
    intronic_floor_fraction = intronic_floor_fraction,
    debris_pct_intronic_prior = debris_pct_intronic_prior,
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

  return(list(
    bounds_empty = bounds_empty,
    bounds_non_empty = bounds_non_empty,
    debris_intronic_floor = debris_intronic_floor
  ))
}

filterNonEmptyPartitionByEmptyIntronicBound <- function(
    df,
    umiThreshold,
    bounds_empty,
    intronic_floor_fraction = 1,
    debris_pct_intronic_prior = 0.25,
    verbose = FALSE) {

  df_high_umi <- df[log10(df$num_transcripts + 1) > umiThreshold, ]

  observed_intronic_floor <- intronic_floor_fraction *
    bounds_empty$intronic_lower_bound

  debris_intronic_floor <- min(
    observed_intronic_floor,
    debris_pct_intronic_prior
  )

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
  local_idx <- getLocalWindowIndex(
    i = i,
    n = length(candidate_results),
    neighbor_window = neighbor_window
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
    localSilhouette = NA_real_,
    stabilityScore = NA_real_,
    selectionScore = NA_real_
) {
  data.frame(
    smoothingMultiplier = smoothingMultiplier,
    umiThreshold = umiThreshold,
    silhouetteScore = silhouetteScore,
    localSilhouette = localSilhouette,
    stabilityScore = stabilityScore,
    selectionScore = selectionScore
  )
}

#' Add a failed smoothing-result row
#'
#' @param resultDF Data frame of per-iteration smoothing results.
#' @param smoothingMultiplier Numeric scalar smoothing multiplier.
#' @param umiThreshold Numeric scalar UMI threshold.
#'
#' @return `resultDF` with one additional row containing missing scores.
#' @noRd
addFailedSmoothingResult <- function(
    resultDF,
    smoothingMultiplier,
    umiThreshold) {

  rbind(
    resultDF,
    makeSmoothingResultRow(
      smoothingMultiplier = smoothingMultiplier,
      umiThreshold = umiThreshold,
      silhouetteScore = NA_real_
    )
  )
}
