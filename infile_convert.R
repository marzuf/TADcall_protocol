startTime <- Sys.time()

topDom_file <- "/media/electron/mnt/etemp/marie/TADcall_yuanlong/25kb/input_caller/chr6/GM12878_chr6_25kb_matrix_pos_zero.txt"
topDom_file <- "//mnt/etemp/marie/TADcall_yuanlong/25kb/input_caller/chr6/GM12878_chr6_25kb_matrix_pos_zero.txt"
topDom_file_test <- "GM12878_chr6_25kb_matrix_pos_zero.txt_test"
catch_file_test <- "GM12878_chr6_25kb_matrix_catch.txt_test"


if(! require(data.table)) {
  install.packages("data.table")
  library(data.table)
} 
if( !require(Matrix)){
  install.packages("Matrix")
  library(Matrix)
} 

TopDom_to_CaTCH <- function(infile, outfile, inSep="\t", outSep="\t", outZeroBased=TRUE) {
  
  cat(paste0("... read: ", infile, "\n"))
  topDom_dt <- fread(infile, sep=inSep, header=FALSE )
  cat(paste0("... done\n"))
  
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
  write.table(out_dt, file=outfile, row.names = FALSE, col.names = FALSE, sep=outSep)
  cat(paste0("... written: ", outfile, "\n"))
  
}
TopDom_to_CaTCH(infile=topDom_file, outfile=catch_file_test)



CaTCH_to_TopDom <- function(infile, outfile, binSize, inSep="\t", outSep="\t", outZeroBased=TRUE, matDim=NULL) {
  
  cat(paste0("... read: ", infile, "\n"))
  catch_dt <- fread(infile, sep=inSep, header=FALSE )
  cat(paste0("... done\n"))
  
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
    binA = seq(0, binSize, length.out=nrow(catch_dt)),
    binB = seq(binSize, binSize, length.out=nrow(catch_dt)),
    stringsAsFactors = FALSE
  )

  topDom_dt <- cbind(threeCols, catch_dt)
    
  stopifnot(ncol(topDom_dt) == nrow(topDom_dt) + 3)
  
  
  
  write.table(topDom_dt, file=outfile, row.names = FALSE, col.names = FALSE, sep=outSep)
  cat(paste0("... written: ", outfile, "\n"))
  
}
CaTCH_to_TopDom(infile=catch_file_test,
                outfile=topDom_file_test,
                binSize=25000,
                matDim = 6843)
  


##########################################################################################
##########################################################################################
##########################################################################################

cat("*** DONE \n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


##########################################################################################
########################################################################################## toy stuff
#####################################################################################