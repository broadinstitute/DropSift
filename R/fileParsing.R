#' Read in DGE matrix
#'
#' Read in a dense DGE matrix in McCarroll lab format (cell barcodes in columns, genes in rows) or
#' read the expected 10x MTX files from a directory.
#'
#' Validates and reorders DGE matrix to match the cell_features order.
#'
#' @param dgeMatrixFile The file to parse.  If a directory, then the expected 10x MTX files are read.
#' @param cell_features A cell features data frame with the cell barcodes to match the DGE matrix.  If supplied,
#' the DGE matrix will be reordered and filtered to match the cell_features barcodes.
#'
#' @return A matrix of the DGE matrix with (optionally) the cell barcodes in the same order as the cell_features
#' @export
readDgeFile<-function (dgeMatrixFile, cell_features=NULL) {
    dgeMatrix=NULL
    if (!is.null(dgeMatrixFile)) {
        dgeMatrix=readDgeMatrixOptionallySparse(dgeMatrixFile)
        #validate the dgeMatrix contains all the cell barcodes in the report
        if (!all(colnames(dgeMatrix) %in% cell_features$cell_barcode)) {
            stop("The cell barcodes in the DGE matrix do not match the cell barcodes in the cell features file.")
        }
        #reorder the dgeMatrix to match the cell_features order
        if (!is.null(cell_features)) {
            dgeMatrix=dgeMatrix[,match(cell_features$cell_barcode, colnames(dgeMatrix))]
            if (any(cell_features$cell_barcode!=colnames(dgeMatrix))) {
                stop ("Cell barcodes in the DGE matrix do not match the cell barcodes in the cell features file.")
            }
        }

    }
    return (dgeMatrix)
}

#' Read in DGE matrix.
#' @inheritParams readDgeFile
#' @noRd
readDgeMatrixOptionallySparse <- function(dgeMatrixFile) {

    #parse the DGE matrix if it is a 10x mtx file.
    dgeMatrix=parse10xMTX(dgeMatrixFile)
    if (!is.null(dgeMatrix))
        return (dgeMatrix)

    #otherwise, read the dense DGE file.
    log_info(paste("Reading expression from dense DGE matrix [", dgeMatrixFile, "]", sep=""))
    # Read the DGE file using DropSeq.utilities
    m <- read_dge_gz(dgeMatrixFile)
    dge <- as(as.matrix(m), "dgCMatrix")
    log_info("Done parsing expression data")
    return(dge)
}

parse10xMTX<-function (dgeMatrixFile) {
    #test to see if the input file is a directory.
    isDirectory <- file.info(dgeMatrixFile)$isdir

    #if this wasn't a directory, then it can't be a 10x mtx file set.
    if (!isDirectory) {
        return (NULL)
    }

    #check if this has the required MTX file names and parse them if so.

    expectedMtxFile=paste(dgeMatrixFile, "/matrix.mtx.gz", sep="")
    expectedFeatureFile=paste(dgeMatrixFile, "/features.tsv.gz", sep="")
    expectedBarcodeFile=paste(dgeMatrixFile, "/barcodes.tsv.gz", sep="")
    mtxExists=file.exists(expectedMtxFile)
    featureExists=file.exists(expectedFeatureFile)
    barcodeExists=file.exists(expectedBarcodeFile)

    #if all the files exist, then read them in.
    if (mtxExists && featureExists && barcodeExists) {
        log_info(paste("Reading expression from MTX file [", expectedMtxFile, "]", sep=""))
        #use the Seurat framework to read the MTX file.
        dge=Seurat::Read10X(data.dir = dirname (expectedMtxFile))
        log_info("Done parsing expression data")
        return(dge)
    } else {
        strErrorMessage=paste("The directory [", dgeMatrixFile, "] does not contain the required files for a 10X mtx file.  Missing files: ",sep="")
        if (!mtxExists) {
            strErrorMessage=paste(strErrorMessage, "matrix.mtx.gz")
        }
        if (!featureExists) {
            strErrorMessage=paste(strErrorMessage, "features.tsv.gz")
        }
        if (!barcodeExists) {
            strErrorMessage=paste(strErrorMessage, "barcodes.tsv.gz")
        }
        stop("The directory does not contain the required files for a 10X mtx file.")
    }

}


#' Read a data.table using data.table fread()
#'
#' If inFile ends with .gz, it is gunzipped in a way that is robust to running out of disk space.
#' Environment variable DROPSEQ_TMPDIR should be set to a directory that has plenty of room and is writable.
#'
#' @param inFile path of file to be read
#' @param comment_regexp If defined, filter input through "egrep -v comment_regexp"
#' @param ... passed through to data.table::fread
#' @return a data.table
#' @import data.table
#' @noRd
fastRead<-function(inFile, comment_regexp=NA, ...) {
    if (!file.exists(inFile)) {
        stop(paste0(inFile, " does not exist"))
    }
    # fread writes process output to /dev/shm if it exists, which is OK so long as it is big enough.
    # If that fails, fall back to function that uses an explicit TMPDIR
    if (length(grep(".gz", inFile))==1) {
        unzip_command = paste("gunzip","-c",inFile,sep=" ")
        if (!is.na(comment_regexp)) {
            unzip_command = paste0(unzip_command, sprintf(" | egrep -v '%s'", comment_regexp))
        }
        # data.table 1.11.8 emits a warning if command string isn't passed via cmd=, but v 1.10.4 (current dotkit version)
        # doesn't recognized cmd=
        if (version_above("data.table", "1.11.8")) {
            return(tryCatch(data.table::fread(cmd=unzip_command,data.table=T, ...), error=function(e) fastReadBigGz(inFile, comment_regexp=comment_regexp, ...)))
        } else {
            return(tryCatch(data.table::fread(unzip_command,data.table=T, ...), error=function(e) fastReadBigGz(inFile, comment_regexp=comment_regexp, ...)))
        }

    } else {
        input = inFile
        if (!is.na(comment_regexp)) {
            input = sprintf("egrep -v '%s' %s", comment_regexp, inFile)
            a=data.table::fread(input, data.table=T, ...)
        } else {
            a=data.table::fread(input, data.table=T, ...)
        }
        return(a)
    }
}

version_above <- function(pkg, than) {
    as.logical(utils::compareVersion(as.character(utils::packageVersion(pkg)), than))
}

# internal function
fastReadBigGz<-function(inFile, ...) {
    # Grid Engine overwrites TMPDIR, so support our own environment variable that won't get clobbered.
    tmpDir = Sys.getenv('DROPSEQ_TMPDIR')
    if (nchar(tmpDir) == 0) {
        tmpDir = Sys.getenv('TMPDIR')
    }
    if (nchar(tmpDir) == 0) {
        tmpDir = tempdir()
    }
    t=tempfile(tmpdir=tmpDir, fileext=".tsv")
    on.exit(unlink(t), add = TRUE)
    cmd=paste("gunzip -c", inFile, ">", t)
    retval = system(cmd,intern=F) != 0
    if (retval != 0) {
        stop(paste0(cmd, " failed with status ", retval))
    }
    return(fastRead(t, ...))
}

#' Read a tabular DGE, with gene rows and columns that typically are cell barcodes
#' @param file to be read, optionally gzipped
#' @param decreasing_order_by_size If true, columns are ordered in decreasing order by colsum
#' @return a data.frame with genes as rownames, no GENE column.
#' @import data.table
#' @noRd
read_dge_gz<-function(file, decreasing_order_by_size=TRUE) {
    dge <- fastRead(file, comment_regexp = '^#', header = TRUE)
    setkey(dge,GENE)
    gene_names = dge$GENE
    GENE <- NULL # Silence R CMD check warning https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
    dge[, GENE:=NULL]
    rownames(dge) = gene_names

    # optionally order cells by expression
    if (decreasing_order_by_size) {
        colsums.full.data=unlist(dge[,lapply(.SD,sum,na.rm=TRUE),.SDcols=colnames(dge)])
        order.data=names(colsums.full.data)[order(colsums.full.data,decreasing=T)]
        full.data.cells=dge[,order.data,with=F]
    } else {
        full.data.cells = dge
    }

    full.data.cells = as.data.frame(full.data.cells)
    rownames(full.data.cells) = gene_names
    return(full.data.cells)
}

#' Load in the cell features file.
#'
#' Loads in the cell features, optionally sorts
#' @param cellFeaturesFile The file to parse.
#' @param  requiredColumns If defined, these column names must be present in the data.  If this check
#' fails the missing required columns are logged and execution is halted.
#' @param verbose If true print message with number of cell features loaded.
#' @return contents of cellFeaturesFile
#' @noRd
readCellFeatures<-function(cellFeaturesFile, requiredColumns=NULL, verbose=TRUE) {
    cell_features <- utils::read.table(cellFeaturesFile, stringsAsFactors = FALSE, header=TRUE )

    if (verbose. && "num_transcripts" %in% colnames(cell_features))
        message( nrow(cell_features), " rows in file (", min(cell_features$num_transcripts), "+ UMIs)")
    return(cell_features)
}




