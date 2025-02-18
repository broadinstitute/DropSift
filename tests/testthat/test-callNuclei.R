test_that("callNuclei", {
    data(svmNucleusCallerInputs)
    set.seed(1)
    svmNucleusCaller = SvmNucleusCaller(
        svmNucleusCallerInputs$cellFeatures,
        svmNucleusCallerInputs$dgeMatrix,
        datasetName = "v3_Bamboo_18d_10X_RNAseq_Optiprep8000_CaH_17k_rxn1",
        useCBRBFeatures = FALSE)
    selectedNuclei = svmNucleusCaller$cell_features[as.logical(svmNucleusCaller$cell_features$is_cell),]$cell_barcode
    expect_equal(length(selectedNuclei), 756)
})
