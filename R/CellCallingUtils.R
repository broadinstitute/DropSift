# MIT License Copyright 2017 Broad Institute Permission is hereby
# granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the 'Software'), to
# deal in the Software without restriction, including without
# limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions: The above copyright notice and this
# permission notice shall be included in all copies or substantial
# portions of the Software.  THE SOFTWARE IS PROVIDED 'AS IS',
# WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#' Get the first pit after the highest peak in the UMI density
#'
#' This function estimates the density of `x`, finds the highest density peak,
#' and returns the first pit to the right of that peak. This is intended to
#' identify the point where the dominant empty-droplet peak ends, rather than
#' a later valley in the high-UMI tail.
#'
#' @param x Numeric vector.
#' @param adjust Numeric scalar passed to `stats::density()`.
#' @param too_close_peak_10X Numeric scalar or `NA`. If supplied, and the next
#'   peak after the highest peak is closer than this distance, the function uses
#'   the pit after that nearby peak. This preserves the existing 10X-specific
#'   behavior.
#' @param density_n Integer scalar. Number of grid points used for density
#'   estimation.
#'
#' @return Numeric scalar x-coordinate of the selected pit, or `NA_real_` if no
#'   valid pit is found.
#' @importFrom pastecs turnpoints
#' @importFrom stats density
#' @noRd
PitAfterHighestPeak <- function(
    x,
    adjust = 2,
    too_close_peak_10X = NA,
    density_n = 200) {

  den <- stats::density(x, adjust = adjust, n = density_n)
  tpts <- pastecs::turnpoints(den$y)

  peaks.wh <- which(tpts$peaks)
  pits.wh <- which(tpts$pits)

  if (length(peaks.wh) == 0 || length(pits.wh) == 0) {
    return(NA_real_)
  }

  peaks.y <- den$y[peaks.wh]
  highpeak.wh1 <- which.max(peaks.y)
  highpeak.index <- peaks.wh[highpeak.wh1]

  candidate_pits <- pits.wh[pits.wh > highpeak.index]

  if (length(candidate_pits) == 0) {
    return(NA_real_)
  }

  nextpit.wh <- candidate_pits[1]

  if (!is.na(too_close_peak_10X) && length(peaks.wh) > highpeak.wh1) {
    nextpeak.index <- peaks.wh[highpeak.wh1 + 1]
    nextpeak.distance <- den$x[nextpeak.index] - den$x[highpeak.index]

    if (nextpeak.distance < too_close_peak_10X) {
      candidate_pits_after_next_peak <- pits.wh[pits.wh > nextpeak.index]

      if (length(candidate_pits_after_next_peak) > 0) {
        nextpit.wh <- candidate_pits_after_next_peak[1]
      }
    }
  }

  den$x[nextpit.wh]
}

#' Find the pit between the two highest peaks
#'
#' This function estimates the density of `x`, identifies the two highest
#' density peaks, and returns the least dense pit between them.
#'
#' @inheritParams PitAfterHighestPeak
#'
#' @return Numeric scalar x-coordinate of the selected pit, or `NA_real_` if the
#'   density does not contain two recognizable peaks or a pit between them.
#' @importFrom pastecs turnpoints
#' @noRd
PitBetweenHighestPeaks <- function(
    x,
    adjust = 2,
    density_n = 200) {

  den <- stats::density(x, adjust = adjust, n = density_n)
  tpts <- pastecs::turnpoints(den$y)

  df <- data.frame(
    index = seq_len(length(tpts$peaks)),
    is_peak = tpts$peaks,
    is_pit = tpts$pits,
    density = den$y[tpts$pos],
    x = den$x[tpts$pos]
  )

  peaks <- df[df$is_peak, ]
  peaks <- peaks[order(peaks$density, decreasing = TRUE), ]

  if (nrow(peaks) < 2) {
    return(NA_real_)
  }

  indexes <- sort(peaks[seq_len(2), ]$index)

  if (any(is.na(indexes))) {
    logger::log_error("Unable to identify peak indexes.")
    return(NA_real_)
  }

  pitBetweenPeaks <- df[
    df$index %in% indexes[1]:indexes[2] & df$is_pit,
  ]

  if (nrow(pitBetweenPeaks) == 0) {
    return(NA_real_)
  }

  pitBetweenPeaks <- pitBetweenPeaks[
    order(pitBetweenPeaks$density, decreasing = FALSE),
  ]

  den$x[pitBetweenPeaks$index[1]]
}

#' Run PitAfterHighestPeak over multiple bandwidths
#'
#' This function runs `PitAfterHighestPeak()` across a grid of density bandwidth
#' adjustment values. The final threshold is selected from a locally stable
#' neighborhood of thresholds across adjacent smoothing values.
#'
#' @inheritParams PitAfterHighestPeak
#' @param adjust Numeric vector of bandwidth adjustment values.
#' @param neighbor_window Integer scalar. Number of neighboring smoothing values
#'   on each side used to compute local threshold stability.
#'
#' @return Numeric scalar UMI threshold.
#' @noRd
PitAfterHighestPeakWithGridSearch <- function(
    x,
    adjust = seq(0.25, 3, 0.25),
    density_n = 200,
    neighbor_window = 2) {

  getThreshold <- function(bw) {
    PitAfterHighestPeak(x, adjust = bw, density_n = density_n)
  }

  thresholds <- vapply(adjust, getThreshold, numeric(1))

  summarizeStablePitThresholds(
    thresholds = thresholds,
    neighbor_window = neighbor_window
  )
}

#' Run PitBetweenHighestPeaks over multiple bandwidths
#'
#' This function runs `PitBetweenHighestPeaks()` across a grid of density
#' bandwidth adjustment values and summarizes the resulting thresholds using
#' their median. This method is used for forced two-cluster solutions.
#'
#' @inheritParams PitBetweenHighestPeaks
#' @param adjust Numeric vector of bandwidth adjustment values.
#'
#' @return Numeric scalar UMI threshold.
#' @noRd
PitBetweenHighestPeaksWithGridSearch <- function(
    x,
    adjust = seq(0.25, 3, 0.25),
    density_n = 200) {

  getThreshold <- function(bw) {
    PitBetweenHighestPeaks(x, adjust = bw, density_n = density_n)
  }

  thresholds <- vapply(adjust, getThreshold, numeric(1))
  summarizePitThresholds(thresholds)
}

#' Summarize pit thresholds with a median
#'
#' @param thresholds Numeric vector of candidate thresholds.
#'
#' @return Numeric scalar threshold, or `NA_real_` if no finite threshold is
#'   available.
#' @noRd
summarizePitThresholds <- function(thresholds) {
  thresholds <- thresholds[is.finite(thresholds)]

  if (length(thresholds) == 0) {
    return(NA_real_)
  }

  stats::median(thresholds, na.rm = TRUE)
}

#' Select a locally stable threshold across smoothing values
#'
#' This function scores each threshold by its local stability across adjacent
#' smoothing values. For each candidate threshold, the local median summarizes
#' the threshold in a neighborhood around that smoothing value, and the local MAD
#' measures how much nearby thresholds vary around that local median.
#'
#' Stability is compared using the local relative MAD, computed as the local MAD
#' divided by the absolute local median. This normalizes the variability before
#' comparing neighborhoods with different threshold scales.
#'
#' @param thresholds Numeric vector of candidate thresholds, ordered by the
#'   smoothing values that generated them.
#' @param neighbor_window Integer scalar. Number of neighboring smoothing values
#'   on each side used to compute local stability.
#'
#' @return Numeric scalar threshold, or `NA_real_` if no finite threshold is
#'   available.
#' @noRd
summarizeStablePitThresholds <- function(
    thresholds,
    neighbor_window = 2) {

  valid_idx <- which(is.finite(thresholds))

  if (length(valid_idx) == 0) {
    return(NA_real_)
  }

  if (length(valid_idx) == 1) {
    return(thresholds[valid_idx])
  }

  local_median <- rep(NA_real_, length(thresholds))
  local_mad <- rep(NA_real_, length(thresholds))
  local_relative_mad <- rep(NA_real_, length(thresholds))

  for (i in valid_idx) {

    idx <- getLocalWindowIndex(
      i = i,
      n = length(thresholds),
      neighbor_window = neighbor_window
    )

    local_values <- thresholds[idx]
    local_values <- local_values[is.finite(local_values)]

    if (length(local_values) == 0) {
      next
    }

    local_median[i] <- stats::median(local_values, na.rm = TRUE)
    local_mad[i] <- stats::median(
      abs(local_values - local_median[i]),
      na.rm = TRUE
    )

    if (local_median[i] != 0) {
      local_relative_mad[i] <- local_mad[i] / abs(local_median[i])
    }
  }

  valid_score_idx <- which(
    is.finite(local_median) & is.finite(local_relative_mad)
  )

  if (length(valid_score_idx) == 0) {
    return(NA_real_)
  }

  local_median <- local_median[valid_score_idx]
  local_relative_mad <- local_relative_mad[valid_score_idx]

  #plot (local_median, local_relative_mad)
  best_idx <- which.min(local_relative_mad)

  local_median[best_idx]

}

#' Get a fixed-width local window around an index
#'
#' This helper returns indices for a local window centered on `i` when possible.
#' Near the beginning or end of the vector, the window is shifted inward rather
#' than truncated. This keeps edge candidates eligible while making their local
#' summaries comparable to candidates in the middle of the vector.
#'
#' @param i Integer scalar. Index to center the window on when possible.
#' @param n Integer scalar. Total number of elements.
#' @param neighbor_window Integer scalar. Number of neighboring elements to use
#'   on each side when the full centered window is available.
#'
#' @return Integer vector of window indices.
#' @noRd
getLocalWindowIndex <- function(i, n, neighbor_window) {
  window_size <- 2 * neighbor_window + 1

  start <- i - neighbor_window
  end <- i + neighbor_window

  if (start < 1) {
    start <- 1
    end <- min(n, window_size)
  }

  if (end > n) {
    end <- n
    start <- max(1, n - window_size + 1)
  }

  seq(start, end)
}
