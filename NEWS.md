# DropSift News

## DropSift (development)
### 🔧 Bug Fixes & Improvements
- `runIntronicSVM()`/`SvmNucleusCaller()` no longer abort when the empty gene
  module score cannot be computed (e.g. when nucleus exemplars are dominated
  by empty droplets and no differentially expressed genes are found). The
  feature is dropped from the SVM and the run continues with a warning.
- Gene module diagnostic plots that could not be computed now render a
  labelled "GENE MODULE SCORE FAILED" placeholder instead of a blank panel.

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
