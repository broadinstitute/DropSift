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

  expect_equal(length(selectedNuclei), 332)
})
