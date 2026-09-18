test_that("safePlotPanel returns the plot on success", {
  p <- DropSift:::safePlotPanel(
    ggplot2::ggplot() +
      ggplot2::theme_void(),
    strTitle = "ok",
    reason = "should not be used"
  )
  expect_s3_class(p, "ggplot")
})

test_that("safePlotPanel returns a placeholder when the expression errors", {
  p <- DropSift:::safePlotPanel(
    stop("boom"),
    strTitle = "broken plot",
    reason = "should not be used"
  )
  expect_s3_class(p, "ggplot")
})

test_that("safePlotPanel returns a placeholder when the expression is NULL", {
  p <- DropSift:::safePlotPanel(
    NULL,
    strTitle = "missing plot",
    reason = "no data available"
  )
  expect_s3_class(p, "ggplot")
})

test_that("getCellSelectionPlotTitle handles zero classified nuclei", {
  df <- data.frame(
    barcode_class = rep(NA_character_, 5),
    num_transcripts = c(10, 20, 30, 40, 50),
    pct_intronic = c(0.1, 0.2, 0.3, 0.4, 0.5)
  )

  strTitle <- DropSift:::getCellSelectionPlotTitle(df, strTitlePrefix = "SVM")

  expect_type(strTitle, "character")
  expect_length(strTitle, 1)
  expect_false(is.na(strTitle))
  expect_true(grepl("0 Nuclei", strTitle))
})

test_that("SvmNucleusCaller degrades gracefully when the nucleus-vs-empty SVM cannot be trained", {
  data(svmNucleusCallerInputs)
  set.seed(1)

  # Force runBinarySVM() to report failure for every binary classifier,
  # exactly as it does in production when there are too few exemplars of
  # one class (see runBinarySVM()'s own class-count guard). This isolates
  # the regression under test (does SvmNucleusCaller() degrade gracefully
  # when the nucleus-vs-empty SVM is unavailable?) from the separate
  # question of what real data triggers that failure, which is reproduced
  # directly against real data in test-runIntronicSVM.R.
  testthat::local_mocked_bindings(
    runBinarySVM = function(...) NULL
  )

  svmNucleusCaller <- expect_no_error(SvmNucleusCaller(
    cellFeatures = svmNucleusCallerInputs$cellFeatures,
    dgeMatrix = svmNucleusCallerInputs$dgeMatrix,
    datasetName = "svm_unavailable_dataset",
    useCBRBFeatures = FALSE,
    forceTwoClusterSolution = FALSE
  ))

  expect_false(svmNucleusCaller$nucleus_svm_available)
  expect_true(length(svmNucleusCaller$degradedWarnings) > 0)
  expect_true(any(grepl(
    "nucleus-vs-empty SVM unavailable",
    svmNucleusCaller$degradedWarnings
  )))
  expect_true(all(is.na(svmNucleusCaller$cell_features$is_cell_prob)))
  expect_false(any(
    svmNucleusCaller$cell_features$barcode_class == "nucleus",
    na.rm = TRUE
  ))

  outFile <- tempfile(fileext = ".pdf")
  grDevices::pdf(outFile)
  expect_no_error(plotSvmNucleusCaller(svmNucleusCaller))
  grDevices::dev.off()
  expect_true(file.exists(outFile))
  # This is the actual observed symptom of the original bug: plotting
  # failed partway through the page and left a truncated/empty PDF.
  expect_gt(file.size(outFile), 0)
})

test_that("getCBRBArgs returns NA when the nucleus-vs-empty SVM was unavailable", {
  data(svmNucleusCallerInputs)
  set.seed(1)

  testthat::local_mocked_bindings(
    runBinarySVM = function(...) NULL
  )

  svmNucleusCaller <- SvmNucleusCaller(
    cellFeatures = svmNucleusCallerInputs$cellFeatures,
    dgeMatrix = svmNucleusCallerInputs$dgeMatrix,
    datasetName = "svm_unavailable_dataset",
    useCBRBFeatures = FALSE,
    forceTwoClusterSolution = FALSE
  )

  cbrbArgs <- expect_no_error(getCBRBArgs(svmNucleusCaller))
  expect_equal(cbrbArgs, list(
    total_droplets_included = NA_integer_,
    expected_cells = NA_integer_
  ))
})

test_that("runIntronicSVM warns instead of erroring when the nucleus-vs-empty SVM cannot be trained", {
  data(svmNucleusCallerInputs)
  set.seed(1)

  testthat::local_mocked_bindings(
    runBinarySVM = function(...) NULL
  )

  tempDir <- tempdir()
  dgeFiles <- writeExampleSvmNucleusCallerInputs(tempDir)
  cellFeaturesFile <- writeExampleCellFeatures(tempDir)
  outPDF <- tempfile(fileext = ".pdf")
  outFeaturesFile <- tempfile(fileext = ".tsv")
  outCbrbFile <- tempfile(fileext = ".tsv")

  expect_warning(
    runIntronicSVM(
      datasetName = "svm_unavailable_dataset",
      cellFeaturesFile = cellFeaturesFile,
      dgeMatrixFile = tempDir,
      useCBRBFeatures = FALSE,
      forceTwoClusterSolution = FALSE,
      outPDF = outPDF,
      outFeaturesFile = outFeaturesFile,
      outCellBenderInitialParameters = outCbrbFile,
      random.seed = 1
    ),
    regexp = "degraded result"
  )

  expect_true(file.exists(outPDF))
  expect_gt(file.size(outPDF), 0)
  expect_true(file.exists(outCbrbFile))
  cbrbResult <- read.delim(outCbrbFile)
  expect_true(is.na(cbrbResult$total_droplets_included))
  expect_true(is.na(cbrbResult$expected_cells))
})
