# DropSift

**DropSift** is a nucleus caller for single-nucleus RNA sequencing data.

This tool distinguishes between nuclei and empty droplets by combining summary metrics of cell barcodes—such as total UMI counts, fraction of intronic reads, and percentage of mitochondrial reads. It also incorporates gene expression to learn signatures characteristic of empty droplets, and can optionally make use of CellBender’s [remove-background](https://github.com/broadinstitute/CellBender) outputs. DropSift uses these features to robustly identify true nuclei across a variety of datasets.

## Installation

To install the latest version directly from GitHub:

```
if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github("broadinstitute/DropSift")
```

## Documentation
A complete guide—including input file formats, example usage, and a detailed walkthrough of the algorithmic steps—is available [here](https://html-preview.github.io/?url=https://github.com/broadinstitute/DropSift/blob/46da7e910692479a5d232e676c382da1eb4f1692/DropSift.html)

## License
This package is released under the MIT License.

