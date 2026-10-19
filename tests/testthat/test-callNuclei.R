test_that("callNuclei", {
  data(svmNucleusCallerInputs)
  set.seed(1)

  svmNucleusCaller <- SvmNucleusCaller(
    svmNucleusCallerInputs$cellFeatures,
    svmNucleusCallerInputs$dgeMatrix,
    datasetName = "v3_Bamboo_18d_10X_RNAseq_Optiprep8000_CaH_17k_rxn1",
    useCBRBFeatures = FALSE
  )

  selectedNuclei <- svmNucleusCaller$cell_features[
    which(svmNucleusCaller$cell_features[["barcode_class"]] == "nucleus"),
  ]$cell_barcode

  expect_equal(length(selectedNuclei), 372)
})

test_that("callNucleiCBRB", {
  data(svmNucleusCallerInputs)
  set.seed(1)

  svmNucleusCaller <- SvmNucleusCaller(
    svmNucleusCallerInputs$cellFeatures,
    svmNucleusCallerInputs$dgeMatrix,
    datasetName = "v3_Bamboo_18d_10X_RNAseq_Optiprep8000_CaH_17k_rxn1",
    useCBRBFeatures = TRUE
  )

  selectedNuclei <- svmNucleusCaller$cell_features[
    which(svmNucleusCaller$cell_features[["barcode_class"]] == "nucleus"),
  ]$cell_barcode

  # The CBRB initialization path always uses a 20 UMI floor for empty
  # exemplar selection, rather than the maxUmisEmpty argument, so this count
  # differs from the non-CBRB test above.
  expect_equal(length(selectedNuclei), 331)
})

test_that("callDebris = FALSE trains only the nucleus-vs-empty SVM", {
  data(svmNucleusCallerInputs)
  set.seed(1)

  expect_no_warning(
    svmNucleusCaller <- SvmNucleusCaller(
      svmNucleusCallerInputs$cellFeatures,
      svmNucleusCallerInputs$dgeMatrix,
      datasetName = "v3_Bamboo_18d_10X_RNAseq_Optiprep8000_CaH_17k_rxn1",
      useCBRBFeatures = FALSE,
      callDebris = FALSE
    )
  )

  cf <- svmNucleusCaller$cell_features
  expect_false(any(cf$barcode_class == "debris", na.rm = TRUE))
  expect_true(all(is.na(cf$is_debris_prob)))
  # low-UMI barcodes are set to NA after the SVMs run.
  idx <- !is.na(cf$is_cell_prob)
  expect_true(any(idx))
  expect_equal(cf$is_cell_prob[idx], cf$p_nucleus_vs_empty[idx])
})
