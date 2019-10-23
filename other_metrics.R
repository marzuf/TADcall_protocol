# 
#' Calculate "Jaccard index" based on concordance of each bin
#'
#' @param file1 Path to file for 1st partition (3-column BED format chromo-start-end; no header) or 3-column dataframe object
#' @param file2 Path to file for 2nd partition (3-column BED format chromo-start-end; no header) or 3-column dataframe object
#' @param binSize Bin size [bp]
#' @param chrSize Chromo size [bp]
#' @param nCpu Nbr of cpu available [use of foreach]


# requirement: foreach, doMC (if nCpu > 1)

get_bin_JaccardIndex <- function(file1, file2, binSize, chrSize=NULL, nCpu = 1){
  
  library(foreach)
  if(nCpu > 1) {
    library(doMC)
    registerDoMC(cores=nCpu)
  }
  
  if(typeof(file1) == "character") {
    if (file.info(file1)$size == 0) {
      set1DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set1DT <- read.delim(file1, header = FALSE, stringsAsFactors = FALSE)
      colnames(set1DT) <- c("chromo", "start", "end")
    }
  } else {
    set1DT <- file1
    stopifnot(ncol(set1DT) == 3)
    colnames(set1DT) <- c("chromo", "start", "end")
  }
  stopifnot(is.numeric(set1DT$start))
  stopifnot(is.numeric(set1DT$end))
  
  if(typeof(file2) == "character") {
    if (file.info(file2)$size == 0) {
      set2DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set2DT <- read.delim(file2, header=FALSE, stringsAsFactors = FALSE)
      colnames(set2DT) <- c("chromo", "start", "end")
    }
  } else {
    set2DT <- file2
    stopifnot(ncol(set2DT) == 3)
    colnames(set2DT) <- c("chromo", "start", "end")
  }
  stopifnot(is.numeric(set2DT$start))
  stopifnot(is.numeric(set2DT$end))

  stopifnot(length(unique(set1DT$chromo)) == 1)
  stopifnot(length(unique(set2DT$chromo)) == 1)
  stopifnot(unique(set2DT$chromo) == unique(set1DT$chromo))
  
  cat(paste0("# of domains in set1\t=\t", nrow(set1DT), "\n"))
  cat(paste0("# of domains in set2\t=\t", nrow(set2DT), "\n"))
  
  
  if(is.null(chrSize)) {
    chrSize <- max(c(set1DT$end, set2DT$end))
  }
  
  # check if the mid-position of the bin is within domain or not
  # "to" should be round up to bin size
  # e.g. if binSize=10 and chrSize=101 -> should go up to 110
  all_bins_mid <- seq(from=binSize/2, to=ceiling(chrSize/binSize)*binSize, by=binSize)
  
  stopifnot(!all_bins_mid %in% set1DT$start)
  stopifnot(!all_bins_mid %in% set1DT$end)
  stopifnot(!all_bins_mid %in% set2DT$start)
  stopifnot(!all_bins_mid %in% set2DT$end)
  pos=all_bins_mid[1]
  matchingVect <- foreach(pos=all_bins_mid, .combine='c') %dopar% {
    as.numeric(any(set1DT$start < pos & set1DT$end > pos) == any(set2DT$start < pos & set2DT$end > pos))
  }
  return(mean(matchingVect))
}
##################################################################################################################################################################################################  
##################################################################################################################################################################################################
##################################################################################################################################################################################################
  
#' Calculate "Jaccard index" for the boundaries 
#'
#' @param file1 Path to file for 1st partition (3-column BED format chromo-start-end; no header) or 3-column dataframe object
#' @param file2 Path to file for 2nd partition (3-column BED format chromo-start-end; no header) or 3-column dataframe object
#' @param tolRad Tolerance radius for boundary matching [bp]
#' @param matchFor Matching for set1 ("set1") or set2 ("set2") only or for both ("all"; default)
#' @param nCpu Nbr of cpu available [use of foreach]

# requirement: GenomicRanges, foreach, doMC (if nCpu > 1)

get_boundaries_JaccardIndex <- function(file1, file2, tolRad, matchFor="all", nCpu = 1){
  
  stopifnot(matchFor %in% c("set1", "set2", "all"))
  
  library(GenomicRanges)
  library(foreach)
  if(nCpu > 1) {
    library(doMC)
    registerDoMC(cores=nCpu)
  }
  
  if(typeof(file1) == "character") {
    if (file.info(file1)$size == 0) {
      set1DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set1DT <- read.delim(file1, header = FALSE, stringsAsFactors = FALSE)
      colnames(set1DT) <- c("chromo", "start", "end")
    }
  } else {
    set1DT <- file1
    stopifnot(ncol(set1DT) == 3)
    colnames(set1DT) <- c("chromo", "start", "end")
  }
  stopifnot(is.numeric(set1DT$start))
  stopifnot(is.numeric(set1DT$end))
  
  if(typeof(file2) == "character") {
    if (file.info(file2)$size == 0) {
      set2DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set2DT <- read.delim(file2, header=FALSE, stringsAsFactors = FALSE)
      colnames(set2DT) <- c("chromo", "start", "end")
    }
  } else {
    set2DT <- file2
    stopifnot(ncol(set2DT) == 3)
    colnames(set2DT) <- c("chromo", "start", "end")
  }
  stopifnot(is.numeric(set2DT$start))
  stopifnot(is.numeric(set2DT$end))
  
  stopifnot(length(unique(set1DT$chromo)) == 1)
  stopifnot(length(unique(set2DT$chromo)) == 1)
  stopifnot(unique(set2DT$chromo) == unique(set1DT$chromo))
  
  mychr <- unique(set1DT$chromo)
  
  cat(paste0("# of domains in set1\t=\t", nrow(set1DT), "\n"))
  cat(paste0("# of domains in set2\t=\t", nrow(set2DT), "\n"))

  # boundaries matching -> treat starts and ends equally:
  # (handle 1-based and 0-based start coordinates)
  if(unique(set1DT$start %%10) == 1) { 
    bd1_DT <- data.frame(chromo=mychr, bdPos=c(set1DT$start-1, set1DT$end))
  } else {
    bd1_DT <- data.frame(chromo=mychr, bdPos=c(set1DT$start, set1DT$end))
  }
  bd1_DT <- unique(bd1_DT)
  bd1_DT <- bd1_DT[order(bd1_DT$bdPos),]
  if(unique(set2DT$start %%10) == 1) { 
    bd2_DT <- data.frame(chromo=mychr, bdPos=c(set2DT$start-1, set2DT$end))
  } else {
    bd2_DT <- data.frame(chromo=mychr, bdPos=c(set2DT$start, set2DT$end))
  }
  bd2_DT <- unique(bd2_DT)
  bd2_DT <- bd2_DT[order(bd2_DT$bdPos),]
  cat(paste0("# of boundaries in set1\t=\t", nrow(bd1_DT), "\n"))
  cat(paste0("# of boundaries in set2\t=\t", nrow(bd2_DT), "\n"))
  
  if(matchFor %in% c("set1", "all")) {
    query1_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bd1_DT$bdPos, end=bd1_DT$bdPos))
    ref2_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bd2_DT$bdPos-tolRad, end=bd2_DT$bdPos+tolRad))
    match1 <- query1_GR %over% ref2_GR
    stopifnot(length(match1) == length(query1_GR))
    stopifnot(length(match1) == nrow(bd1_DT))
    if(matchFor == "set1") return(mean(match1))
  }
  if(matchFor %in% c("set2", "all")) {
    query2_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bd2_DT$bdPos, end=bd2_DT$bdPos))
    ref1_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bd1_DT$bdPos-tolRad, end=bd1_DT$bdPos+tolRad)) 
    match2 <- query2_GR %over% ref1_GR
    stopifnot(length(match2) == length(query2_GR))
    stopifnot(length(match2) == nrow(bd2_DT))
    if(matchFor == "set2") return(mean(match2))
  }
  # if come to here -> matchFor=="all"
  return(mean(c(match1, match2)))
}

##################################################################################################################################################################################################  
##################################################################################################################################################################################################
##################################################################################################################################################################################################

#' Calculate ratio of matching TAD
#'
#' @param file1 Path to file for 1st partition (3-column BED format chromo-start-end; no header) or 3-column dataframe object
#' @param file2 Path to file for 2nd partition (3-column BED format chromo-start-end; no header) or 3-column dataframe object
#' @param coverMatchRatioThresh Ratio to be covered to be considered matching [0-1]
#' @param matchFor Matching for set1 ("set1") or set2 ("set2") only or for both ("all"; default)
#' @param nCpu Nbr of cpu available [use of foreach]
#' 
#' WARNING: coverMatchRatioThresh=0 -> will not return 1 -> will return ratio of TADs that have a match, without considering a threshold of bp matching !

# requirement: GenomicRanges, foreach, doMC (if nCpu > 1)

get_ratioMatchingTADs <- function(file1, file2, coverMatchRatioThresh, matchFor="all", nCpu = 1){
  
  stopifnot(matchFor %in% c("set1", "set2", "all"))
  
  stopifnot(coverMatchRatioThresh >= 0 & coverMatchRatioThresh <= 1)
  
  library(GenomicRanges)
  library(foreach)
  if(nCpu > 1) {
    library(doMC)
    registerDoMC(cores=nCpu)
  }
  
  if(typeof(file1) == "character") {
    if (file.info(file1)$size == 0) {
      set1DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set1DT <- read.delim(file1, header = FALSE, stringsAsFactors = FALSE)
      colnames(set1DT) <- c("chromo", "start", "end")
    }
  } else {
    set1DT <- file1
    stopifnot(ncol(set1DT) == 3)
    colnames(set1DT) <- c("chromo", "start", "end")
  }
  stopifnot(is.numeric(set1DT$start))
  stopifnot(is.numeric(set1DT$end))
  
  if(typeof(file2) == "character") {
    if (file.info(file2)$size == 0) {
      set2DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set2DT <- read.delim(file2, header=FALSE, stringsAsFactors = FALSE)
      colnames(set2DT) <- c("chromo", "start", "end")
    }
  } else {
    set2DT <- file2
    stopifnot(ncol(set2DT) == 3)
    colnames(set2DT) <- c("chromo", "start", "end")
  }
  stopifnot(is.numeric(set2DT$start))
  stopifnot(is.numeric(set2DT$end))
  
  stopifnot(length(unique(set1DT$chromo)) == 1)
  stopifnot(length(unique(set2DT$chromo)) == 1)
  stopifnot(unique(set2DT$chromo) == unique(set1DT$chromo))
  
  mychr <- unique(set1DT$chromo)
  
  cat(paste0("# of domains in set1\t=\t", nrow(set1DT), "\n"))
  cat(paste0("# of domains in set2\t=\t", nrow(set2DT), "\n"))
  
  # (handle 1-based and 0-based start coordinates) # for TAD size calc.
  if(unique(set1DT$start %%10) == 1) set1DT$start <- set1DT$start-1
  if(unique(set2DT$start %%10) == 1) set2DT$start <- set2DT$start-1

  
  set1DT$region <- paste0("set1_", mychr, "_TAD", 1:nrow(set1DT))
  set1_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=set1DT$start, end=set1DT$end, names=set1DT$region))
  
  set2DT$region <- paste0("set2_", mychr, "_TAD", 1:nrow(set2DT))
  set2_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=set2DT$start, end=set2DT$end, names=set2DT$region))
  
  
  if(matchFor %in% c("set1", "all")) {
    ref_GR <- set1_GR
    query_GR <- set2_GR
    ref_tad_size_set1 <- setNames(c(set1DT$end-set1DT$start), c(set1DT$region))
    # determine which features from the query overlap which features in the subject
    overlap_GR <- findOverlaps(query=query_GR, subject=ref_GR)
    if(length(overlap_GR) == 0) {  
      set1_match <- 0
    } else {
      IDoverlap_hits_all <- findOverlaps(query=query_GR,
                                         subject=query_GR)
      IDoverlap_hits <- IDoverlap_hits_all[queryHits(IDoverlap_hits_all) != subjectHits(IDoverlap_hits_all)]
      
      refID <- names(ref_GR[subjectHits(overlap_GR)])
      queryID <- names(query_GR[queryHits(overlap_GR)])
      
      stopifnot(refID %in% set1DT$region)
      stopifnot(refID %in% names(ref_tad_size_set1))
      
      set1_overlapDT_bp <- data.frame(
        refID = refID,
        queryID = queryID,
        overlapBp = width(pintersect(ref_GR[refID], query_GR[queryID])),
        overlapBpRatio = width(pintersect(ref_GR[refID], query_GR[queryID]))/ref_tad_size_set1[refID],
        stringsAsFactors = FALSE)
      # retain only the matching that pass the threshold
      set1_overlapDT_bp <- set1_overlapDT_bp[set1_overlapDT_bp$overlapBpRatio >= coverMatchRatioThresh,]
      set1_match <- set1DT$region %in% set1_overlapDT_bp$refID
      stopifnot(length(set1_match) == nrow(set1DT))
    }
    if(matchFor == "set1") return(mean(set1_match))
  }
  
  if(matchFor %in% c("set2", "all")) {
    ref_GR <- set2_GR
    query_GR <- set1_GR
    ref_tad_size_set2 <- setNames(c(set2DT$end-set2DT$start), c(set2DT$region))
    # determine which features from the query overlap which features in the subject
    overlap_GR <- findOverlaps(query=query_GR, subject=ref_GR)
    if(length(overlap_GR) == 0) {  
      set2_match <- 0
    } else {
      IDoverlap_hits_all <- findOverlaps(query=query_GR,
                                         subject=query_GR)
      IDoverlap_hits <- IDoverlap_hits_all[queryHits(IDoverlap_hits_all) != subjectHits(IDoverlap_hits_all)]
      
      refID <- names(ref_GR[subjectHits(overlap_GR)])
      queryID <- names(query_GR[queryHits(overlap_GR)])
      
      stopifnot(refID %in% set2DT$region)
      stopifnot(refID %in% names(ref_tad_size_set2))
      
      set2_overlapDT_bp <- data.frame(
        refID = refID,
        queryID = queryID,
        overlapBp = width(pintersect(ref_GR[refID], query_GR[queryID])),
        overlapBpRatio = width(pintersect(ref_GR[refID], query_GR[queryID]))/ref_tad_size_set2[refID],
        stringsAsFactors = FALSE)
      # retain only the matching that pass the threshold
      set2_overlapDT_bp <- set2_overlapDT_bp[set2_overlapDT_bp$overlapBpRatio >= coverMatchRatioThresh,]
      set2_match <- set2DT$region %in% set2_overlapDT_bp$refID
      stopifnot(length(set2_match) == nrow(set2DT))
    }
    if(matchFor == "set2") return(mean(set2_match))
  }
  # if come to here -> matchFor=="all"
  return(mean(c(set1_match, set2_match)))
  
}



######################################################################## trial toy stuff
# 
# tolRad = 10*2
# 
# bddt1 = data.frame(chr="chr6",
#                  bdPos = c(0, 30, 50,100,140,180))
# 
# bddt2  = data.frame(chr="chr6",
#                    bdPos = c(10, 80, 140, 190, 220))
# 
# 
# query1_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bddt1$bdPos, end=bddt1$bdPos))
# ref1_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bddt1$bdPos-tolRad, end=bddt1$bdPos+tolRad))
# 
# query2_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bddt2$bdPos, end=bddt2$bdPos))
# ref2_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bddt2$bdPos-tolRad, end=bddt2$bdPos+tolRad))
# 
# 
# query1_GR %over% ref1_GR
# query2_GR %over% ref2_GR
# 
# 
# query1_GR %over% ref2_GR
# query2_GR %over% ref1_GR








                   