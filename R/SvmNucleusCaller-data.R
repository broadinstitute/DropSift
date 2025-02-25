# MIT License Copyright 2024 Broad Institute Permission is hereby
# granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the 'Software'), to
# deal in the Software without restriction, including without
# limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions: The above copyright notice and this
# permission notice shall be included in all copies or substantial
# portions of the Software.  THE SOFTWARE IS PROVIDED 'AS IS',
# WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#' svmNucleusCallerInputs
#'
#' Example cell features and digital gene expression (DGE) matrix for testing
#' the SVM nucleus caller. Data is from a single-nucleus RNA sequencing
#' (snRNA-seq) experiment:
#' **2023-08-25_v3_Bamboo_18d_10X_RNAseq_Optiprep8000_CaH_17k_rxn1**.
#'
#' @docType data
#'
#' @usage data(svmNucleusCallerInputs)
#'
#' @format A list with two elements:
#' \describe{
#'   \item{cellFeatures}{A data frame containing cell-level features.}
#'   \item{dgeMatrix}{A sparse matrix of digital gene expression (DGE)
#'   data.}
#' }
#'
#' @keywords datasets
#'
#' @details
#' The data used in this vignette is sourced from a single sequencing reaction
#' from the [BICAN](https://assets.nemoarchive.org/col-npe62xg) project.
#'
#' **Citation:**
#' Kim K, Macaisa L, Drouin S, Rayan N, Mello C, Yavari N, Yoo O, Wysoker A,
#' Shakir K, Nemesh J, Kashin S, Goldman M, Genovese G, Fritch H, Hogan M,
#' Flowers K, Finn E, Vanderburg C, Ichihara K, Macosko E, McCarroll S (2024).
#' *Population variability: Single nucleus RNAseq from postnatal human brain*
#' Dataset. Available from:
#' [NEMO Archive](https://assets.nemoarchive.org/col-npe62xg).
#'
#' @examples
#' data(svmNucleusCallerInputs)
#' svmNucleusCaller <- DropSift::SvmNucleusCaller(
#'   svmNucleusCallerInputs$cellFeatures,
#'   svmNucleusCallerInputs$dgeMatrix,
#'   datasetName = "v3_Bamboo_18d_10X_RNAseq_Optiprep8000_CaH_17k_rxn1",
#'   useCBRBFeatures = FALSE
#' )
#'
"svmNucleusCallerInputs"
