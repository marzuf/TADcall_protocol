startTime <- Sys.time()

options(scipen=100)

topDom_file <- "/media/electron/mnt/etemp/marie/TADcall_yuanlong/25kb/input_caller/chr6/GM12878_chr6_25kb_matrix_pos_zero.txt"
topDom_file <- "//mnt/etemp/marie/TADcall_yuanlong/25kb/input_caller/chr6/GM12878_chr6_25kb_matrix_pos_zero.txt"
topDom_file_test <- "GM12878_chr6_25kb_matrix_pos_zero.txt_test"
catch_file_test <- "GM12878_chr6_25kb_matrix_catch.txt_test"
arrowhead_file_test <- "GM12878_chr6_25kb_matrix.pre"


if(! require(data.table)) {
  install.packages("data.table")
  library(data.table)
}
if( !require(Matrix)){
  install.packages("Matrix")
  library(Matrix)
}

run_TopDom_to_CaTCH <- FALSE
run_CaTCH_to_TopDom <- FALSE
run_TopDom_to_arrowhead <- TRUE

##############################################################################################
############################################################################################## => TopDom -> CaTCH
##############################################################################################

TopDom_to_CaTCH <- function(infile, outfile, inSep="\t", outSep="\t", outZeroBased=FALSE) {
  
  cat(paste0("... read: ", infile, "\n"))
  topDom_dt <- fread(infile, sep=inSep, header=FALSE )
  cat(paste0("... done\n"))
  topDom_dt <- data.frame(topDom_dt)
  chromo <- unique(as.character(topDom_dt[,1]))
  stopifnot(length(chromo) == 1)
  
  stopifnot(ncol(topDom_dt) == nrow(topDom_dt) + 3)
  topDom_dt <- topDom_dt[,-c(1:3)]
  
  stopifnot(ncol(topDom_dt) == nrow(topDom_dt))
  
  topDom_mat <- as(as.matrix(topDom_dt), "sparseMatrix") 
  
  i <- topDom_mat@i
  j <- rep(seq_along(diff(topDom_mat@p)), diff(topDom_mat@p))
  count_values <- topDom_mat@x
  stopifnot(length(i) == length(j))
  stopifnot(length(i) == length(count_values))
  
  if(outZeroBased) {
    i <- i-1
    j <- j-1
  }

  # CaTCH OUTUT FORMAT:  
  # col1 = chromosome
  # col2 = bin of the start region (genomic coordinate divided by binsize)
  # col3 = bin of the end region (genomic coordinate divided by binsize)
  # col4 = Hi-C counts
  
  out_dt <- data.frame(
    chromo = chromo,
    binStart = i,
    binEnd = j,
    count = count_values,
    stringsAsFactors = FALSE
  ) 
  write.table(out_dt, file=outfile, row.names = FALSE, col.names = FALSE, sep=outSep, quote=FALSE)
  cat(paste0("... written: ", outfile, "\n"))
  
}
if(run_TopDom_to_CaTCH) TopDom_to_CaTCH(infile=topDom_file, outfile=catch_file_test)

##############################################################################################
############################################################################################## => CaTCH -> TopDom 
##############################################################################################


CaTCH_to_TopDom <- function(infile, outfile, binSize, inSep="\t", outSep="\t", inZeroBased=FALSE, matDim=NULL) {
  
  cat(paste0("... read: ", infile, "\n"))
  catch_dt <- fread(infile, sep=inSep, header=FALSE )
  cat(paste0("... done\n"))
  catch_dt <- data.frame(catch_dt)
  stopifnot(ncol(catch_dt) == 4)
  colnames(catch_dt) <- c("chromo", "binA", "binB", "count")

  
  chromo <- unique(as.character(catch_dt[,1]))
  stopifnot(length(chromo) == 1)

  i <- catch_dt$binA + as.numeric(inZeroBased) # 1-based for SparseMatrix
  j <- catch_dt$binB + as.numeric(inZeroBased)
  stopifnot(i > 0)
  stopifnot(j > 0)
  count_values <- catch_dt$count
  
  if(is.null(matDim)) {
    matDim <- max(c(i,j))
  }  else {
    stopifnot(matDim >= c(i,j))
  }
  catch_mat <- sparseMatrix(i = i, j = j, x = count_values, dims=c(matDim,matDim))  # default: the index vectors i and/or j are 1-based

  catch_dt <- as.data.frame(as.matrix(catch_mat))
  stopifnot(nrow(catch_dt) == ncol(catch_dt))
  
  threeCols <- data.frame(
    chromo = rep(chromo, nrow(catch_dt)),
    binA = seq(0, by=binSize, length.out=nrow(catch_dt)),
    binB = seq(binSize, by=binSize, length.out=nrow(catch_dt)),
    stringsAsFactors = FALSE
  )

  topDom_dt <- cbind(threeCols, catch_dt)
    
  stopifnot(ncol(topDom_dt) == nrow(topDom_dt) + 3)
  
  write.table(topDom_dt, file=outfile, row.names = FALSE, col.names = FALSE, sep=outSep)
  cat(paste0("... written: ", outfile, "\n"))
  
}
if(run_CaTCH_to_TopDom) 
  CaTCH_to_TopDom(infile=catch_file_test,
                outfile=topDom_file_test,
                binSize=25000,
                matDim = 6843)
  
##############################################################################################
############################################################################################## => TopDom -> arrowhead (pre) 
##############################################################################################


TopDom_to_arrowhead <- function(infile, outfile, inSep="\t",  chrSize=NULL, sizefile=NULL) {

  cat(paste0("... read: ", infile, "\n"))
  topDom_dt <- fread(infile, sep=inSep, header=FALSE )
  topDom_dt <- data.frame(topDom_dt)
  cat(paste0("... done\n"))
  
  chromo <- unique(as.character(topDom_dt[,1]))
  stopifnot(length(chromo) == 1)
  
  stopifnot(is.numeric(topDom_dt[,2]))
  stopifnot(is.numeric(topDom_dt[,3]))
  stopifnot(topDom_dt[,2] < topDom_dt[,3])
  if(is.null(chrSize)) {
    chrSize <- max(topDom_dt[,3])
    stopifnot(chrSize == topDom_dt[nrow(topDom_dt), 3])
  }
  
  binSize <- unique(topDom_dt[,3]-topDom_dt[,2])
  stopifnot(length(binSize) == 1)
  
  stopifnot(ncol(topDom_dt) == nrow(topDom_dt) + 3)
  topDom_dt <- topDom_dt[,-c(1:3)]
  
  if(is.null(sizefile)) {
    sizefile <- file.path(dirname(outfile), paste0(chromo, ".size"))
  }
  sizedt <- data.frame(chr=chromo, size=chrSize, stringsAsFactors = FALSE)
  write.table(sizedt, file=sizefile, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE) # !!! should be tab-separated
  cat(paste0("... written: ", sizefile, "\n"))
  
  topDom_mat <- Matrix(as.matrix(topDom_dt), sparse=TRUE)
  matrixDF <- data.frame(summary(topDom_mat))
  stopifnot(ncol(matrixDF) == 3)
  # the strand can be both 0 everywhere, but not the fragment
  # nRows <-nrow(matrixDF)
  fragVec1 <- rep(0,  nrow(matrixDF))
  fragVec2 <- rep(1,  nrow(matrixDF))
  posVec1 <- (matrixDF[,1] - 1 ) * binSize # 0-based
  posVec2 <- (matrixDF[,2] - 1 ) * binSize # 0-based
  countVec <- matrixDF[,3]
  # chrVec <- rep(sub("chr", "", chromo), nRows)
  # chrVec <- rep(chromo, nRows)
  # preDT <- data.frame(str1=fragVec1, chr1=chrVec, pos1=posVec1, frag1=fragVec1, 
  #                     str1=fragVec1, chr1=chrVec, pos1=posVec2, frag1=fragVec2, counts=countVec)
  # <frag1> should be 0 and <frag2> should be 1.
  # <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2>
  preDT <- data.frame(str1=0, chr1=chromo, pos1=pmin(posVec1, posVec2), frag1=0, 
                      str1=0, chr1=chromo, pos1=pmax(posVec1, posVec2), frag1=1, counts=countVec)
  write.table(preDT, file=outfile, row.names=FALSE, col.names=FALSE, sep=" ", quote=FALSE) # !!! should be white space-separated
  cat(paste0("... written: ", outfile, "\n"))
  
}
if(run_TopDom_to_arrowhead)
  TopDom_to_arrowhead(infile=topDom_file,
                outfile=arrowhead_file_test)

##########################################################################################
##########################################################################################
##########################################################################################

cat("*** DONE \n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


##########################################################################################
########################################################################################## toy stuff
##########################################################################################
# catch_dt <- data.frame(chromo="chr1", binA=c(1,1,1,2,3,3), binB=c(1,2,3,2,2,3), count=c(10,2,1,12,7,9))
# 
# countM <- matrix(c(1,3,4,2,3,8,7,5,4,7,10,9,2,5,9,11), byrow=T, nrow=4)
# sparseMatrix <- Matrix(as.matrix(countM), sparse=T)
# matrixDF <- data.frame(summary(sparseMatrix))
# 
# countM <- matrix(c(1,0,4,2,0,8,0,5,4,0,10,9,2,5,9,11), byrow=T, nrow=4)
# sparseMatrix <- Matrix(as.matrix(countM), sparse=T)
# matrixDF <- data.frame(summary(sparseMatrix))
# compare CaTCH format

