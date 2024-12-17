# MIT License
# 
# Copyright 2024 Broad Institute
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#' svmNucleusCallerInputs
#'
#' cell features and digital gene expression matrix for testing the SVM nucleus caller.
#' Data from a snRNA-seq experiment 2023-08-25_v3_Bamboo_18d_10X_RNAseq_Optiprep8000_CaH_17k_rxn1_dev_aug07_2024
#'
#' @docType data
#'
#' @usage data(svmNucleusCallerInputs)
#'
#' @format A list with 2 elements:
#' cellFeatures: a data frame of cell features
#' dgeMatrix: a sparse matrix of digital gene expression
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(svmNucleusCallerInputs)
#' svmNucleusCaller = Dropseq.nucleiselection::SvmNucleusCaller(svmNucleusCallerInputs$cellFeatures,
#' svmNucleusCallerInputs$dgeMatrix, datasetName = "v3_Bamboo_18d_10X_RNAseq_Optiprep8000_CaH_17k_rxn1",
#' useCBRBFeatures = FALSE)
#'
"svmNucleusCallerInputs"
