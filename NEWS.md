# DropSift News

## DropSift (development)
### 🔧 Bug Fixes & Improvements
- `runIntronicSVM()`/`SvmNucleusCaller()` no longer abort when the
  nucleus-vs-empty SVM itself cannot be trained (e.g. too few nucleus
  exemplars). The run now completes with an R `warning()` describing the
  degradation instead of stopping before any output is written; the
  resulting cell features have `barcode_class`/`is_cell_prob` set to `NA`
  for every barcode, and `getCBRBArgs()` returns `NA` CellBender parameters
  instead of erroring.
- `runIntronicSVM()`/`SvmNucleusCaller()` no longer abort when the empty gene
  module score cannot be computed (e.g. when nucleus exemplars are dominated
  by empty droplets and no differentially expressed genes are found). The
  feature is dropped from the SVM and the run continues with a warning.
- Individual SVM cell-selection plot panels that fail to build or render (for
  example because no barcode was classified as nucleus) now render a
  labelled "PLOT UNAVAILABLE" placeholder in that panel's position instead of
  aborting the whole page; the PDF device is now always closed via
  `on.exit()`, and the SVM training-selection panel is emitted as a
  last-resort fallback if the full plot layout fails unexpectedly. This
  fixes a bug where a degraded run could produce a 0-byte/truncated PDF.
- Gene module diagnostic plots that could not be computed now render a
  labelled "GENE MODULE SCORE FAILED" placeholder instead of a blank panel.
  These placeholders now also respect `useCBRBFeatures`, so a degraded run
  with `useCBRBFeatures = FALSE` no longer shows placeholders for CellBender
  plots that were never requested.
- `runIntronicSVM()`/`SvmNucleusCaller()`/`findTrainingDataBounds()` gain a
  new `twoClusterFallbackRatio` parameter (default `2`). When the default
  exemplar-selection solution's silhouette score is weak, the two-cluster
  solution is now computed automatically and adopted in place of the default
  when its silhouette score is at least `twoClusterFallbackRatio` times
  better, correcting cases where the default solution's nucleus exemplars
  were diluted by a near-ambient density mode. This replaces the need to
  manually decide when to pass `forceTwoClusterSolution = TRUE`. Set
  `twoClusterFallbackRatio = NULL` to disable and always use the default
  solution.

## DropSift 1.0.0 (Bioconductor Release)
### 🆕 Initial Release
- First public release of `DropSift` on Bioconductor.
- Implements **SVM-based nuclei selection** from single-nucleus RNA-seq (snRNA-seq) data.
- Supports **input from 10x Genomics, Optimus H5AD, and dense DGE formats**.
- **Automated feature extraction**, including:
  - UMI counts
  - % intronic reads
  - Mitochondrial content
  - Empty gene module scores
- **Flexible classifier initialization**, handling:
  - With/without CellBender remove-background.
  - High ambient RNA contamination cases.
- **Efficient sparse matrix processing** with `Matrix` and `data.table`.
- **Comprehensive visualization tools**:
  - Quality control plots.
  - Feature distributions for classifier training.
  - Selection probability visualization.
- **Test dataset (`svmNucleusCallerInputs`)** included for reproducible examples.

### 🔧 Bug Fixes & Improvements
- Improved **input validation & error handling**.
- Ensured **Bioconductor compliance**, including:
  - Fully runnable examples.
  - Properly formatted documentation.
