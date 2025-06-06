test_that("runIntronicSVM", {
    # Create a temporary directory for test data
    temp_dir <- tempdir()

    # Write example data to files
    dge_files <- writeExampleSvmNucleusCallerInputs(temp_dir)
    cell_features_file <- writeExampleCellFeatures(temp_dir)

    # Define output file paths
    out_pdf_file <- file.path(temp_dir, "output.pdf")
    out_features_file <- file.path(temp_dir, "classified_features.tsv")
    out_cbrb_args_file <- file.path(temp_dir, "cbrb_args.tsv")

    #if set.seed isn't set, the result can vary slightly.
    set.seed(1)
    # Run the function
    runIntronicSVM(
       datasetName = "example_dataset",
       cellFeaturesFile = cell_features_file,
       dgeMatrixFile = temp_dir,  # directory with compressed 10x files
       useCBRBFeatures = FALSE,
       forceTwoClusterSolution = FALSE,
       outPDF = out_pdf_file,
       outFeaturesFile = out_features_file,
       outCellBenderInitialParameters=out_cbrb_args_file
    )

    # Check that output files were generated
    file.exists(out_pdf_file)    # Should return TRUE
    file.exists(out_features_file)  # Should return TRUE
    file.exists(out_cbrb_args_file)  # Should return TRUE

    result<- read.table(out_features_file, sep = "\t",
        header = TRUE, stringsAsFactors = FALSE)

    selectedNuclei <-  result[which(result$is_cell==T),]$cell_barcode
    expect_equal(length(selectedNuclei), 376)

    cbrb_args_result<- read.table(out_cbrb_args_file, sep = "\t",
        header = TRUE, stringsAsFactors = FALSE)
    expect_equal(cbrb_args_result$total_droplets_included, 2195)
    expect_equal(cbrb_args_result$expected_cells, 376)
})
