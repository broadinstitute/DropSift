test_that("parse h5", {
  h5ad_file <- system.file("extdata", "adata_example.h5ad.gz",
    package = "DropSift"
  )
  r <- parseH5ad(h5ad_file)
  dge <- r$dge
  cell_features <- r$cell_features

  # Check dimensions of the parsed expression matrix
  d <- dim(dge)
  # The example DGE is 5 cells and 5 genes.
  expect_equal(d[1], 5)
  expect_equal(d[2], 5)

  # This cell feature set is the raw outputs from Optimus which are not yet
  # read to be used by DropSift.  See parseOptimusH5ad.
  expect_equal(nrow(cell_features), 5)
  expect_equal(ncol(cell_features), 43)
})

test_that("parse Optimus h5ad", {
  h5ad_file <- system.file("extdata", "adata_example.h5ad.gz",
    package = "DropSift"
  )
  r <- parseOptimusH5ad(h5ad_file, min_transcripts = 0)
  dge <- r$dge
  cell_features <- r$cell_features

  # Check dimensions of the parsed expression matrix
  d <- dim(dge)
  # The example DGE is 5 cells and 5 genes.
  expect_equal(d[1], 5)
  expect_equal(d[2], 5)

  # This cell feature set is the fully processed data for DropSift.
  expect_equal(nrow(cell_features), 5)
  expect_equal(ncol(cell_features), 5)
})
