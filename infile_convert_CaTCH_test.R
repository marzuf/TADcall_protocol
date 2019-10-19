startTime <- Sys.time()

topDom_file <- "/media/electron/mnt/etemp/marie/TADcall_yuanlong/25kb/input_caller/chr6/GM12878_chr6_25kb_matrix_pos_zero.txt"
topDom_file <- "//mnt/etemp/marie/TADcall_yuanlong/25kb/input_caller/chr6/GM12878_chr6_25kb_matrix_pos_zero.txt"
catch_file_test <- "GM12878_chr6_25kb_matrix_catch.txt_test"
catch_file_test_chromo <- "GM12878_chr6_25kb_matrix_catch.txt_test_chr"
catch_file_test_symm <- "GM12878_chr6_25kb_matrix_catch.txt_test_symm"
catch_file_test_symm_sorted <- "GM12878_chr6_25kb_matrix_catch.txt_test_symm_sorted"
catch_file_test_zerobased <- "GM12878_chr6_25kb_matrix_catch.txt_test_zerobased"
catch_file_test_fullsymm <- "GM12878_chr6_25kb_matrix_catch.txt_test_fullsymm"

if(! require(data.table)) {
  install.packages("data.table")
  library(data.table)
}
if( !require(Matrix)){
  install.packages("Matrix")
  library(Matrix)
}

run_TopDom_to_CaTCH <- TRUE

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
    chromo = 1,
    binStart = i,
    binEnd = j,
    count = count_values,
    stringsAsFactors = FALSE
  ) 
  write.table(out_dt, file=outfile, row.names = FALSE, col.names = FALSE, sep=outSep, quote=FALSE)
  cat(paste0("... written: ", outfile, "\n"))
  
}
# if(run_TopDom_to_CaTCH) TopDom_to_CaTCH(infile=topDom_file, outfile=catch_file_test)

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
# if(run_TopDom_to_CaTCH) TopDom_to_CaTCH(infile=topDom_file, outfile=catch_file_test_chromo)

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
  topDom_mat <- forceSymmetric(topDom_mat)
  
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
# if(run_TopDom_to_CaTCH) TopDom_to_CaTCH(infile=topDom_file, outfile=catch_file_test_symm)

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
  
  # CaTCH OUTPUT FORMAT:  
  # col1 = chromosome
  # col2 = bin of the start region (genomic coordinate divided by binsize)
  # col3 = bin of the end region (genomic coordinate divided by binsize)
  # col4 = Hi-C counts
  
  out_dt <- data.frame(
    chromo = 1,
    binStart = i,
    binEnd = j,
    count = count_values,
    stringsAsFactors = FALSE
  ) 
  write.table(out_dt, file=outfile, row.names = FALSE, col.names = FALSE, sep=outSep, quote=FALSE)
  cat(paste0("... written: ", outfile, "\n"))
  
}
# if(run_TopDom_to_CaTCH) TopDom_to_CaTCH(infile=topDom_file, outfile=catch_file_test_zerobased, outZeroBased=TRUE)
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
  topDom_mat <- forceSymmetric(topDom_mat)
  
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
  out_dt <- out_dt[order(out_dt$binStart, out_dt$binEnd),]
  
  write.table(out_dt, file=outfile, row.names = FALSE, col.names = FALSE, sep=outSep, quote=FALSE)
  cat(paste0("... written: ", outfile, "\n"))
  
}
# if(run_TopDom_to_CaTCH) TopDom_to_CaTCH(infile=topDom_file, outfile=catch_file_test_symm_sorted)
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
  topDom_mat <- forceSymmetric(topDom_mat)
  topDom_mat <- as(as.matrix(topDom_mat), "sparseMatrix") 
  
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
  out_dt <- out_dt[order(out_dt$binStart, out_dt$binEnd),]
  
  write.table(out_dt, file=outfile, row.names = FALSE, col.names = FALSE, sep=outSep, quote=FALSE)
  cat(paste0("... written: ", outfile, "\n"))
  
}
if(run_TopDom_to_CaTCH) TopDom_to_CaTCH(infile=topDom_file, outfile=catch_file_test_fullsymm)



##########################################################################################
##########################################################################################
##########################################################################################

cat("*** DONE \n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


##########################################################################################
########################################################################################## toy stuff
##########################################################################################



catch_out_test <- get(load("CaTCH_out_test.Rdata"))
catch_out_test_chr <- get(load("CaTCH_out_test_chr.Rdata"))
catch_out_test_symm <- get(load("CaTCH_out_test_symm.Rdata")) 
catch_out_test_symm_sorted <- get(load("CaTCH_out_test_symm_sorted.Rdata"))
catch_out_test_zerobased <- get(load("CaTCH_out_test_zerobased.Rdata"))
catch_out_test_fullsymm <- get(load("CaTCH_out_test_fullsymm.Rdata"))


dt <- catch_out_test[["clusters"]][,-1]
dt_chr <- catch_out_test_chr[["clusters"]][,-1]
dt_symm <- catch_out_test_symm[["clusters"]][,-1]
dt_symm_sorted <- catch_out_test_symm_sorted[["clusters"]][,-1]
dt_zerobased <- catch_out_test_zerobased[["clusters"]][,-1]
dt_fullsymm <- catch_out_test_fullsymm[["clusters"]][,-1]


all.equal(dt,dt_chr)  # TRUE
all.equal(dt,dt_symm) # FALSE
all.equal(dt,dt_symm_sorted) # FALSE
all.equal(dt,dt_zerobased) # FALSE (but same insulation)
all.equal(dt,dt_fullsymm) # TRUE

all.equal(dt_symm,dt_symm_sorted) # TRUE

all.equal(dt_symm,dt_fullsymm) # FALSE



