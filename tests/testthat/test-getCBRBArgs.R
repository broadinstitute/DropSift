test_that("CBRB argument estimation", {
  data(svmNucleusCallerInputs)
  set.seed(1)
  svmNucleusCaller <- SvmNucleusCaller(
    svmNucleusCallerInputs$cellFeatures,
    svmNucleusCallerInputs$dgeMatrix,
    datasetName = "v3_Bamboo_18d_10X_RNAseq_Optiprep8000_CaH_17k_rxn1",
    useCBRBFeatures = FALSE,
    forceTwoClusterSolution = FALSE
  )
  cbrbArgs <- getCBRBArgs(svmNucleusCaller)
  expect_equal(cbrbArgs, list(
    total_droplets_included = 2195,
    expected_cells = 376
  ))
})
