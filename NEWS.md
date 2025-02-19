# DropSift News

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

## Future Roadmap
- Testing on non-human organisms

