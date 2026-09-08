test_that("makeEmptyGeneModuleResult marks the result invalid with a reason", {
  cell_features_labeled <- data.frame(
    training_label_class = c("nucleus", "empty", "nucleus", "empty")
  )

  result <- DropSift:::makeEmptyGeneModuleResult(
    cell_features_labeled = cell_features_labeled,
    module_score_name = "empty_gene_module_score",
    reason = "no differentially expressed genes found"
  )

  expect_false(result$valid)
  expect_equal(result$reason, "no differentially expressed genes found")
  expect_true(all(is.na(result$score)))
  expect_equal(length(result$score), nrow(cell_features_labeled))

  expectedNames <- c(
    "empty_gene_module_score_training_data",
    "empty_gene_module_score",
    "frac_contamination",
    "empty_gene_module_score_vs_contam"
  )
  expect_true(all(expectedNames %in% names(result$plots)))
  for (plotName in expectedNames) {
    expect_s3_class(result$plots[[plotName]], "ggplot")
  }
})

test_that("computeSvmGeneModuleScore returns valid=FALSE when there are no DE genes", {
  set.seed(1)
  numGenes <- 20
  numCells <- 40

  # A constant relative expression profile across all cells means every
  # gene has an identical normalized value in the nucleus and empty groups,
  # so no gene passes the log-fold-change threshold used by de_wilcox().
  base_prob <- rep(1 / numGenes, numGenes)
  cell_totals <- rep(1000, numCells)
  dgeMatrix <- outer(base_prob, cell_totals)
  dgeMatrix <- round(dgeMatrix)
  rownames(dgeMatrix) <- paste0("GENE-", seq_len(numGenes))
  colnames(dgeMatrix) <- paste0("CELL-", seq_len(numCells))

  cell_features_labeled <- data.frame(
    training_label_class = rep(c("nucleus", "empty"), each = numCells / 2),
    row.names = colnames(dgeMatrix)
  )

  result <- computeSvmGeneModuleScore(
    cell_features_labeled = cell_features_labeled,
    dgeMatrix = dgeMatrix,
    numGenes = 100,
    useCellBenderFeatures = FALSE,
    negative_class = "empty",
    min_nucleus_exemplars = 5,
    min_negative_exemplars = 5,
    min_pseudobulk_observations = 1,
    verbose = FALSE
  )

  expect_false(result$valid)
  expect_true(nzchar(result$reason))
  expect_true(all(is.na(result$score)))
})

test_that("plotGeneModuleScoresByExemplarClass returns a placeholder when scores are unavailable", {
  cell_features_labeled <- data.frame(
    empty_gene_module_score = c(NA_real_, NA_real_, NA_real_),
    debris_gene_module_score = c(0.1, -0.2, 0.3),
    training_label_class = c("empty", "nucleus", "debris")
  )

  p <- plotGeneModuleScoresByExemplarClass(
    cell_features_labeled,
    strTitle = "Gene module scores by exemplar class"
  )

  expect_s3_class(p, "ggplot")
})

test_that("SvmNucleusCaller degrades gracefully when the empty gene module score is unavailable", {
  data(svmNucleusCallerInputs)
  set.seed(1)

  # Force computeSvmGeneModuleScore() to report failure for every gene
  # module, exactly as it does in production when there are no
  # differentially expressed genes. This isolates the regression under test
  # (does SvmNucleusCaller() degrade gracefully?) from the separate question
  # of what real data triggers that failure, which is covered directly by
  # the "computeSvmGeneModuleScore returns valid=FALSE" test above.
  testthat::local_mocked_bindings(
    computeSvmGeneModuleScore = function(cell_features_labeled, ...,
                                         negative_class = "empty") {
      makeEmptyGeneModuleResult(
        cell_features_labeled = cell_features_labeled,
        module_score_name = paste0(negative_class, "_gene_module_score"),
        reason = "no differentially expressed genes found"
      )
    }
  )

  svmNucleusCaller <- SvmNucleusCaller(
    cellFeatures = svmNucleusCallerInputs$cellFeatures,
    dgeMatrix = svmNucleusCallerInputs$dgeMatrix,
    datasetName = "constant_expression_dataset",
    useCBRBFeatures = FALSE,
    forceTwoClusterSolution = FALSE
  )

  expect_false("empty_gene_module_score" %in% svmNucleusCaller$features)
  expect_true(length(svmNucleusCaller$degradedWarnings) > 0)
  expect_true(any(grepl(
    "empty gene module score unavailable",
    svmNucleusCaller$degradedWarnings
  )))
  expect_true(nrow(svmNucleusCaller$cell_features) > 0)

  outFile <- tempfile(fileext = ".pdf")
  grDevices::pdf(outFile)
  expect_no_error(plotSvmNucleusCaller(svmNucleusCaller))
  grDevices::dev.off()
  expect_true(file.exists(outFile))
})
