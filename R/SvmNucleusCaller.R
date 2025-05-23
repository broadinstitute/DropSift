# MIT License Copyright 2024 Broad Institute Permission is hereby
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

# NucleiCaller S3 class definition

# constructor, only used internally
new_SvmNucleusCaller <- function(results, cellProbabilityThreshold,
    maxUmisEmpty, forceTwoClusterSolution, useCBRBInitialization) {
    stopifnot(is.list(results))
    results$useCBRBInitialization <- useCBRBInitialization
    results$cellProbabilityThreshold <- cellProbabilityThreshold
    results$maxUmisEmpty <- maxUmisEmpty
    results$forceTwoClusterSolution <- forceTwoClusterSolution
    return(structure(results, class = c("SvmNucleusCaller", "list")))
}

emptyGeneModuleScoreColName <- "empty_gene_module_score"
requiredNonSvmColNames <- c("cell_barcode")
cbrbNonSvmColNames <- c("num_retained_transcripts")
contaminationColName <- "frac_contamination"

#' Create an SvmNucleusCaller object
#'
#' Train an SVM to distinguish between nuclei and empty droplets.  It may also
#' be used to estimate CellBender remove-background arguments.
#'
#' @param cellFeatures A data frame with columns: cell_barcode (unless
#'   dgeMatrix==NULL), num_reads, frac_contamination (unless !useCBRBFeatures)
#'   and any addtional columns that are specified in featureColumns.
#' @param dgeMatrix A matrix of gene expression values, with rows corresponding
#'   to genes and columns corresponding to cells. The column names must match
#'   the cell_barcode column in cellFeatures.  This may be a sparse or dense
#'   matrix.  Set to NULL to suppress use of gene expression features.
#' @param cellProbabilityThreshold The probability threshold for selecting
#'   cells.  This value scales from 0-1, with 1 being extremely
#'   confident/stringent.  Set to null to use the SVM defaults (>0.5 = nuclei).
#'   This modifies plotting outputs and the output features file classification,
#'   but does not filter any cells from the output.
#' @param maxUmisEmpty All cell barcodes with fewer than this number of UMIs are
#'   not considered for any aspects of nuclei selection.  It is assumed these
#'   cell barcodes have fewer UMIs than the ambient RNA.
#' @param featureColumns A list of features to use for cell selection.  If NULL,
#'   featureColumns are determined based on other parameters.  See
#'   [configureFeatureColumns()].
#' @param forceTwoClusterSolution When true, the initialization of the SVM will
#'   attempt to find a solution with two clusters. In cases where an experiment
#'   is overloaded with nuclei, this may correct the initial set of nuclei and
#'   empty droplets selected. This argument is specific to the density-style
#'   selection of exemplars that is only applicable when useCBRBFeatures is
#'   false.
#' @param useCBRBFeatures When true, the cell bender feature frac_contamination
#'   is used for cell selection. When false, these features are not used.  This
#'   modifies the featureColumns argument.
#' @param useCBRBInitialization When true, the cellbender feature
#'   frac_contamination is used to select exemplar nuclei and empty droplets.
#'   This option can be false and useCBRBFeatures to force DropSift to use the
#'   non-CBRB initialization but still include CBRB in the model features.
#' @param datasetName A string to identify the dataset in plots.
#' @return An SvmNucleusCaller object
#' @seealso [configureFeatureColumns()]
#' @import logger
#' @export
#' @examples
#' # Load example dataset
#' data(svmNucleusCallerInputs)
#' set.seed(1)  # Set seed for reproducibility
#'
#' # Create an SvmNucleusCaller object
#' svmNucleusCaller <- SvmNucleusCaller(
#'     cellFeatures = svmNucleusCallerInputs$cellFeatures,
#'     dgeMatrix = svmNucleusCallerInputs$dgeMatrix,
#'     datasetName = "example_dataset",
#'     useCBRBFeatures = FALSE,
#'     useCBRBInitialization=FALSE,
#'     forceTwoClusterSolution = FALSE
#' )
#' # View summary
#' print(svmNucleusCaller)

SvmNucleusCaller <- function(cellFeatures, dgeMatrix,
    cellProbabilityThreshold = NULL, maxUmisEmpty = 50, featureColumns = NULL,
    forceTwoClusterSolution = FALSE, useCBRBFeatures = TRUE,
    useCBRBInitialization=useCBRBFeatures, datasetName = "") {

    if (useCBRBInitialization==TRUE & useCBRBFeatures==FALSE)
        stop ("Can't use CBRB for initialization of CBCB without features")

    if (useCBRBFeatures==TRUE & useCBRBInitialization==FALSE) {
        msg="Using non-CBRB initialization and CBRB features in model training"
        log_info(msg)
    }

    stopifnot(is.data.frame(cellFeatures))
    stopifnot(is.null(dgeMatrix) || is.matrix(dgeMatrix) || is(dgeMatrix,
        "sparseMatrix"))
    stopifnot(is.null(cellProbabilityThreshold) ||
        is.numeric(cellProbabilityThreshold))

    stopifnot(is.numeric(maxUmisEmpty))
    stopifnot(is.logical(forceTwoClusterSolution))
    stopifnot(is.character(datasetName))
    stopifnot(is.logical(useCBRBFeatures))

    if (!is.null(dgeMatrix)) {
        log_info(sprintf("%d cell barcodes overlap between cellFeatures and dgeMatrix",
            length(intersect(cellFeatures$cell_barcode,
                colnames(dgeMatrix)))))
        stopifnot(nrow(dgeMatrix) > 0)
        stopifnot(ncol(dgeMatrix) > 0)
        stopifnot(all(colnames(dgeMatrix) %in% cellFeatures$cell_barcode))
    }

    # optionally enhance the cell features data frame with additional
    # columns if they are missing and can be computed.
    cellFeatures <- enhanceCelLFeatures(cellFeatures,
        dgeMatrix, useCBRBFeatures)

    featureColumns <- configureFeatureColumns(featureColumns, useCBRBFeatures,
        dgeMatrix)
    stopifnot(is.character(featureColumns))
    validateCellFeatures(cellFeatures, featureColumns)
    if (contaminationColName %in% featureColumns &&
        !all(cbrbNonSvmColNames %in% colnames(cellFeatures))) {

        warning(sprintf("If using CBRB features and column %s is not in ",
            "cellFeatures, plotting will fail.",
            paste(cbrbNonSvmColNames, collapse = ", ")))
    }
    results <- callByIntronicSVM(dataset_name = datasetName,
        cell_features = cellFeatures, dgeMatrix = dgeMatrix,
        cellProbabilityThreshold = cellProbabilityThreshold,
        max_umis_empty = maxUmisEmpty, features = featureColumns,
        useCBRBInitialization=useCBRBInitialization,
        forceTwoClusterSolution = forceTwoClusterSolution)

    results$cell_features <- data.frame(
        cell_barcode = rownames(results$cell_features),
        results$cell_features
    )
    return(new_SvmNucleusCaller(results, cellProbabilityThreshold, maxUmisEmpty,
        forceTwoClusterSolution, useCBRBInitialization))
}

#' Optionally add extra cell features to the cell_features data frame if they
#' are missing and can be computed from the input data.
#' @noRd
enhanceCelLFeatures <- function(cellFeatures, dgeMatrix, useCBRBFeatures) {
    # optionally add the num_transcripts column to the cell_features
    # data frame if missing.
    if (!"num_transcripts" %in% colnames(cellFeatures)) {
        log_info("num_transcripts column not found in cell_features.",
            "Computing from DGE matrix.")
        cellFeatures$num_transcripts <- colSums(dgeMatrix)
    }

    # optionally approximate the num_retained_transcripts if missing
    if (!"num_retained_transcripts" %in% colnames(cellFeatures)
        && useCBRBFeatures) {
        log_info("frac_contamination found but num_retained_transcripts column",
            "not found in cell_features.  Approximating from num_transcripts ",
            "and frac_contamination.")
        cellFeatures$num_retained_transcripts <- round(
            (1 - cellFeatures$frac_contamination) * cellFeatures$num_transcripts
        )

    }
    return(cellFeatures)
}

validateCellFeatures <- function(cellFeatures, features) {
    requiredColumns <- setdiff(features, emptyGeneModuleScoreColName)
    requiredColumns <- c(requiredColumns, requiredNonSvmColNames)
    if (!is.null(requiredColumns)) {
        missingColumns <-
            requiredColumns[!requiredColumns %in% colnames(cellFeatures)]
        if (length(missingColumns) > 0) {
            log_error("The cell features file is missing the following",
            "required columns: ", paste(missingColumns, collapse = ", "))
            stop()
        }
    }
}

DefaultFeatureColumns <- c("num_transcripts", "pct_intronic", "pct_mt")
defaultFeatureColumnsStringRep <- paste0("c(\"", paste(DefaultFeatureColumns,
    collapse = "\", \""), "\")")

#' Send a standard set of plots (3 pages) to the current graphics device
#'
#' @param svmNucleusCaller an object of class SvmNucleusCaller
#' @return the input object, invisibly
#' @export
#' @examples
#' # Load example dataset
#' data(svmNucleusCallerInputs)
#' set.seed(1)  # Set seed for reproducibility
#'
#' # Create an SvmNucleusCaller object
#' svmNucleusCaller <- SvmNucleusCaller(
#'     cellFeatures = svmNucleusCallerInputs$cellFeatures,
#'     dgeMatrix = svmNucleusCallerInputs$dgeMatrix,
#'     datasetName = "example_dataset",
#'     useCBRBFeatures = FALSE,
#'     forceTwoClusterSolution = FALSE
#' )
#'
#' # View summary
#' print(svmNucleusCaller)
#'
#' # Plot the results, directing the outputs to a temp file.
#' outFile=tempfile(fileext=".pdf")
#' pdf(outFile)
#' plotSvmNucleusCaller(svmNucleusCaller)
#' dev.off()

plotSvmNucleusCaller <- function(svmNucleusCaller) {
    UseMethod("plotSvmNucleusCaller", svmNucleusCaller)
}

#' Configure featureColumns argument for `SvmNucleusCaller`
#'
#' Typically, the user will not need to call this function directly.
#'
#' If `featureColumns` is not `NULL`, the user has explicitly specified the
#' features to use. If `featureColumns` is `NULL`, the feature columns will be
#' set to the default values stored in `defaultFeatureColumnsStringRep`.
#' If `useCBRBFeatures` is `TRUE`, the feature column defined by
#' `contaminationColName` will be added. If `dgeMatrix` is not `NULL`,
#' `emptyGeneModuleScoreColName` will be added.
#'
#' @inheritParams SvmNucleusCaller
#' @keywords internal
#' @return The list of features to use to train the SVM
configureFeatureColumns <- function(featureColumns,
    useCBRBFeatures,
    dgeMatrix) {

    if (is.null(featureColumns)) {
        featureColumns <- DefaultFeatureColumns
        if (useCBRBFeatures) {
            featureColumns <- c(featureColumns, contaminationColName)
        }
        if (!is.null(dgeMatrix)) {
            featureColumns <- c(featureColumns, emptyGeneModuleScoreColName)
        }
    }
    return(featureColumns)
}

#' @param svmNucleusCaller an object of class SvmNucleusCaller
#' @rdname plotSvmNucleusCaller
#' @export
plotSvmNucleusCaller.SvmNucleusCaller <- function(svmNucleusCaller) {
    plots <- svmNucleusCaller$plots
    geneModulePlots <- svmNucleusCaller$geneModulePlots
    featurePlot <- svmNucleusCaller$featurePlot
    datasetName <- svmNucleusCaller$dataset_name
    if (contaminationColName %in% svmNucleusCaller$features) {
        arrangeSVMCellSelectionPlots(plots, geneModulePlots = geneModulePlots,
            featurePlot = featurePlot, datasetName,
            outPDF = NULL, useOpenPDF = TRUE)
    } else {
        arrangeSVMCellSelectionPlotsNoCBRB(plots,
            geneModulePlots = geneModulePlots,
            featurePlot = featurePlot, datasetName,
            outPDF = NULL, useOpenPDF = TRUE)
    }
    invisible(svmNucleusCaller)
}

#' @noRd
#' @export
plotSvmNucleusCaller.default <- function(svmNucleusCaller) {
    stop("plotSvmNucleusCaller not implemented for this class")
}

#' Print the SvmNucleusCaller object summary
#' @param x an object of class SvmNucleusCaller
#' @param ... Other arguments that are ignored
#' @method print SvmNucleusCaller
#' @return Invisibly returns the input `SvmNucleusCaller` object.
#' @export
#' @examples
#' # Load example dataset
#' data(svmNucleusCallerInputs)
#' set.seed(1)  # Set seed for reproducibility
#'
#' # Create an SvmNucleusCaller object
#' svmNucleusCaller <- SvmNucleusCaller(
#'     cellFeatures = svmNucleusCallerInputs$cellFeatures,
#'     dgeMatrix = svmNucleusCallerInputs$dgeMatrix,
#'     datasetName = "example_dataset",
#'     useCBRBFeatures = FALSE,
#'     forceTwoClusterSolution = FALSE
#' )
#'
#' # View summary
#' print(svmNucleusCaller)
print.SvmNucleusCaller <- function(x, ...) {
    svmNucleusCaller <- x
    if (is.null(svmNucleusCaller$cellProbabilityThreshold)) {
        cellProbabilityThreshold <- "NULL"
    } else {
        cellProbabilityThreshold <- svmNucleusCaller$cellProbabilityThreshold
    }
    cat("SvmNucleusCaller object\n")
    cat("Dataset: ", svmNucleusCaller$dataset_name, "\n")
    cat("cellProbabilityThreshold: ", cellProbabilityThreshold, "\n")
    cat("maxUmisEmpty: ", svmNucleusCaller$maxUmisEmpty, "\n")
    cat("forceTwoClusterSolution: ", svmNucleusCaller$forceTwoClusterSolution,
        "\n")
    cat("featureColumns: ", paste(svmNucleusCaller$features, collapse = ", "),
        "\n")
    cat("Input cell features: ", nrow(svmNucleusCaller$cell_features),
        " cells\n")
    cat("Number of cells selected: ",
        length(which(as.logical(svmNucleusCaller$cell_features$is_cell))),"\n")
    invisible(svmNucleusCaller)
}

#' Estimate CellBender remove-background arguments from an SvmNucleusCaller
#' object
#'
#' After running the SVM with useCBRBFeatures=FALSE, get an estimate for
#' --total-droplets-included and --expected-cells CellBender remove-background
#' arguments.
#'
#' @param svmNucleusCaller an object of class
#'   SvmNucleusCaller(useCBRBFeatures=FALSE)
#' @return a list with two elements: total_droplets_included and expected_cells
#' @examples
#' #Load up the example dataset and run the SVM.
#' #This mirrors the included unit test.
#' data(svmNucleusCallerInputs)  # Load example dataset
#' set.seed(1)  # Set seed for reproducibility
#' svmNucleusCaller <- SvmNucleusCaller(
#'     svmNucleusCallerInputs$cellFeatures,
#'     svmNucleusCallerInputs$dgeMatrix,
#'     datasetName = "example_dataset",
#'     useCBRBFeatures = FALSE,
#'     forceTwoClusterSolution = FALSE
#' )
#'
#' cbrbArgs <- getCBRBArgs(svmNucleusCaller)
#' print(cbrbArgs)
#'
#' @export
getCBRBArgs <- function(svmNucleusCaller) {
    UseMethod("getCBRBArgs", svmNucleusCaller)
}

#' @param svmNucleusCaller an object of class
#'   SvmNucleusCaller(useCBRBFeatures=FALSE)
#' @export
#' @rdname getCBRBArgs
getCBRBArgs.SvmNucleusCaller <- function(svmNucleusCaller) {
    if (contaminationColName %in% svmNucleusCaller$features) {
        stop("getCBRBArgs should only be used when useCBRBFeatures is false.")
    }
    df <- svmNucleusCaller$cell_features
    threshold_total_droplets <- round(mean(df[df$training_label_is_cell ==
        FALSE, ]$num_transcripts, na.rm = TRUE))
    total_droplets_included <-
        length(which(df$num_transcripts > threshold_total_droplets))

    expected_cells <- length(which(df$is_cell == TRUE))
    return(list(
        total_droplets_included = total_droplets_included,
        expected_cells = expected_cells)
    )

}

#' @export
getCBRBArgs.default <- function(svmNucleusCaller) {
    stop("getCBRBArgs not implemented for this class")
}
