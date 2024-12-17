# MIT License
#
# Copyright 2017 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



#' Get the pit after the highest peak in the UMI density
#'
#' @param x The data
#' @param adjust how much to adjust the bandwith.
#' @param too_close_peak_10X A fix for 10x genomics data.
#' @param density_n Sets the number of bins for density estimation.
#' @return The highest peak
#' @importFrom pastecs turnpoints
#' @importFrom stats density
#' @noRd
PitAfterHighestPeak = function(x, adjust=2, too_close_peak_10X=NA, density_n=200) {
    den = stats::density(x, adjust=adjust, n=density_n )
    #this can fail if N is not set in some circumstances.
    # and setting N slightly changes the results because the standard N is 512 (200 is more course grained)
    #the warning message is pastecs::turnpoints(den$y) : value out of range in 'gammafn'
    #den = density( log10(x), adjust=adjust)
    tpts = pastecs::turnpoints( den$y )
    peaks.wh = which( tpts$peaks )
    pits.wh  = which( tpts$pits  )
    peaks.y = den$y[ peaks.wh ]
    highpeak.wh1 = which.max(peaks.y)

    # This result can be NA when there are bandwidth problems, or when the peaks are inverted.
    nextpit.wh = pits.wh[ highpeak.wh1 ]
    # Don't enable 10X functionality if any of the following are true:
    # * no next peak after the highest peak
    # * no pit after the next peak
    # * next peak is too far from highest peak
    if (!is.na(too_close_peak_10X) && length(peaks.y) > highpeak.wh1 && length(pits.wh) > highpeak.wh1 &&
        den$x[peaks.wh[highpeak.wh1 + 1]] - den$x[peaks.wh[highpeak.wh1]] < too_close_peak_10X) {
        nextpit.wh = pits.wh[ highpeak.wh1+1 ]
    }
    nextpit.x = den$x[nextpit.wh]
    return(nextpit.x)
}

#' Reworked pit after highest peak.
#'
#' Instead of finding the pit after the highest peak, find the pit between the two highest peaks.
#' This finds the two peaks with the highest density, then the pit with the lowest density between them.
#'
#' #Unfortunately, this doesn't work well when there is noise in the data below the empty peak.  Maybe there
#' were some way to identify nuclei vs empty and not get fooled, this might be useful.  Needs work
#' at best, so don't use!
#'
#' @inheritParams PitAfterHighestPeak
#' @return The x value of the pit between the two highest peaks.
#' @importFrom pastecs turnpoints
#' @noRd
PitBetweenHighestPeaks<-function (x, adjust=2, density_n=200) {
    den = density(x, adjust=adjust, n=density_n )
    #this can fail if N is not set in some circumstances.
    # and setting N slightly changes the results because the standard N is 512 (200 is more course grained)
    #the warning message is pastecs::turnpoints(den$y) : value out of range in 'gammafn'
    #den = density( log10(x), adjust=adjust)
    tpts = pastecs::turnpoints( den$y )
    peaks.wh = which( tpts$peaks )
    pits.wh  = which( tpts$pits  )
    peaks.y = den$y[ peaks.wh ]

    #set up the results as a dataframe which is easier to look at.
    #there are some edge cases where turningpoints doesn't return all inputs (ie: less than density_n inputs)
    #so we need to make sure that the indexes are correct.
    df=data.frame(index=1:length(tpts$peaks), is_peak=tpts$peaks, is_pit=tpts$pits, density=den$y[tpts$pos], x=den$x[tpts$pos])

    #select the peaks, order them by density, and select the two largest.
    peaks=df[df$is_peak,]
    peaks=peaks[order(peaks$density, decreasing=TRUE),]
    #if there aren't two recognizable peaks, return NA.
    #hopefully a different density will pick up a reasonable value.
    if (dim(peaks)[1] < 2) {
        return(NA)
    }

    indexes=sort(peaks[1:2,]$index)
    if (any(is.na(indexes))) {
        log_error("This is bad")
    }

    #select the potential pits between these two peak indexes.
    pitBetwenPeaks=df[df$index %in% indexes[1]:indexes[2] & df$is_pit==T,]

    #select the least dense pit between the two peaks.
    pitBetwenPeaks=pitBetwenPeaks[order(pitBetwenPeaks$density, decreasing=TRUE),]
    nextpit.wh = pitBetwenPeaks$index[1]
    nextpit.x = den$x[nextpit.wh]
    return(nextpit.x)
}

#' Runs PitAfterHighestPeakSmarter, but searches multiple bandwidths for non-NA results.
#'
#' PitAfterHighestPeakSmarter works for most but not all data sets.  Changing the level of smoothing
#' of the bandwidth can help.
#' Most of the outputs are very similar for non NA results, so can take the median result.
#' @inheritParams PitAfterHighestPeak
#' @return The median of the non-NA results.
#' @noRd
PitAfterHighestPeakWithGridSearch<-function (x, adjust=seq(0.25, 3, 0.25), density_n=200) {
    f<-function (bw) {
        PitAfterHighestPeak(x, adjust=bw, density_n=density_n)
    }
    r=sapply(adjust, f)
    return (stats::median(r, na.rm=T))
}

#' Runs PitBetweenHighestPeaks, but searches multiple bandwidths for non-NA results.
#'
#' @inheritParams PitBetweenHighestPeaks
#' @return The median of the non-NA results.
#' @noRd
PitBetweenHighestPeaksWithGridSearch<-function (x, adjust=seq(0.25, 3, 0.25), density_n=200) {
    f<-function (bw) {
        PitBetweenHighestPeaks(x, adjust=bw, density_n=density_n)
    }
    r=sapply(adjust, f)
    return (stats::median(r, na.rm=T))
}

