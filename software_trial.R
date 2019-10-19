options(scipen=100)

###########################################################################################################################
# TopDom - package version
###########################################################################################################################
topDom_file <- "/mnt/etemp/marie/TADcall_yuanlong/25kb/input_caller/chr6/GM12878_chr6_25kb_matrix_pos_zero.txt"
stopifnot(file.exists(topDom_file))
topDom_outFile <- file.path("out_TopDom", "topDom_results")
topDom_TADs_outFile <- file.path("out_TopDom", "topDom_final_domains.txt")
dir.create(dirname(topDom_outFile), recursive = TRUE)
dir.create(dirname(topDom_TADs_outFile), recursive = TRUE)

if(!require(TopDom)) {
  if( !require(remotes)) {
    install.packages("remotes")
    library(remotes)
  }
  remotes::install_github("HenrikBengtsson/TopDom")
  library(topDom)
}
##### 1) run TopDom  
topDom_out <- TopDom(topDom_file, window.size=5, outFile=topDom_outFile)
  # > head(topDom_out[["domain"]])
  # chr from.id from.coord to.id to.coord    tag   size
  # 1  chr6       1          0     6   150000    gap 150000
  # 9  chr6       7     150000    16   400000 domain 250000
  # > head(topDom_out[["bed"]])                                                                                                                                                                                 
  # chrom chromStart chromEnd   name
  # 1   chr6          0   150000    gap
  # 9   chr6     150000   400000 domain
  # > head(topDom_out[["binSignal"]])                                                                                                                                                                           
  # id  chr from.coord to.coord local.ext  mean.cf pvalue
  # 1  1 chr6          0    25000      -0.5   0.0000      1
  # 2  2 chr6      25000    50000      -0.5   0.0000      1

##### 2) prepare output
topDom_tads_dt <- topDom_out[["bed"]]
topDom_tads_dt <- topDom_tads_dt[as.character(topDom_tads_dt$name) == "domain",]
out_dt <- topDom_tads_dt[,c("chrom", "chromStart","chromEnd")]
# to 1-based start positions:
out_dt$chromStart <- out_dt$chromStart + 1 
# ensure ordering
out_dt <- out_dt[order(out_dt$chromStart, out_dt$chromEnd),]
write.table(out_dt, file=topDom_TADs_outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
cat(paste0("... written: ", topDom_TADs_outFile, "\n"))

###########################################################################################################################
# TopDom - original script version
###########################################################################################################################

# original site not longer available
# script still available: https://github.com/HenrikBengtsson/TopDom/blob/0.0.2/R/TopDom.R

source("TopDom.R")
# => then use of TopDom() function similar as above

###########################################################################################################################
# CaTCH
###########################################################################################################################
source("infile_convert.R")
topDom_file <- "/mnt/etemp/marie/TADcall_yuanlong/25kb/input_caller/chr6/GM12878_chr6_25kb_matrix_pos_zero.txt"
catch_file <- "GM12878_chr6_25kb_matrix_list.txt"
TopDom_to_CaTCH(infile=topDom_file, outfile=catch_file, outZeroBased=FALSE)

# installation instructions: https://github.com/zhanyinx/CaTCH_R/
# once installed
stopifnot(require(CaTCH))
stopifnot(file.exists(catch_file))
catch_TADs_outFile <- file.path("out_CaTCH", "CaTCH_final_domains.txt")
dir.create(dirname(catch_TADs_outFile), recursive = TRUE)

chromo <- "chr6"
ri_thresh <- 0.65
bin_size <- 25*1000

##### 1) run CaTCH
CaTCH_out <- domain.call(catch_file)

##### 2) prepare output
CaTCH_dt <- CaTCH_out[["clusters"]]
# > head(CaTCH_dt)
# chromosome    RI start end insulation
# 1          1 0.001     1   2        NaN
# 2          1 0.001     3   6  0.3074472
# 3          1 0.001     7   8  0.1807085
# > head(CaTCH_ndt)
# chromosome    RI ndomains
# 1          1 0.000     3421
# 2          1 0.001     3383

CaTCH_ndt <- CaTCH_out[["ncluster"]]
all_ri <- unique(CaTCH_dt$RI)
ref_ri <- all_ri[which.min(abs(all_ri - ri_thresh))]
# ref_ri = 0.7 -> 306 domains
# ref_ri = 0.65 -> 389 domains
# ref_ri = 0.6 -> 526 domains
cat(paste0("... selected RI thresh\t=\t", ri_thresh,"\n"))
cat(paste0("... best matching RI thresh\t=\t", ref_ri,"\n"))
domainsDT <- CaTCH_dt[CaTCH_dt$RI == ref_ri, c(1, 3, 4)]
colnames(domainsDT) <- c("chr", "start", "end")

nDomains <- CaTCH_ndt$ndomains[CaTCH_ndt$RI == ref_ri]
stopifnot(length(nDomains) == 1)
stopifnot(nDomains == nrow(domainsDT))

domainsDT$start <- (domainsDT$start - 1) * bin_size + 1 
domainsDT$end <- domainsDT$end * bin_size 
# convert
# chromosome start end
# 1     1   2
# 1     3   8
# 1     9  23

# to
# chromosome   start     end
# 1       1   50000
# 1   50001  200000
# 1  200001  575000
# ensure ordering
domainsDT <- domainsDT[order(domainsDT$start, domainsDT$end),]
write.table(domainsDT, file=catch_TADs_outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
cat(paste0("... written: ", catch_TADs_outFile), "\n")


# multi-threads (if > 1 chromo) domain.call.parallel(inputs,ncpu=parallel::detectCores()-1L)

###########################################################################################################################
# arrowhead
###########################################################################################################################

chromo <- "chr6"
bin_size <- 25*1000
juicer_pre_file = "GM12878_chr6_25kb_matrix.pre" # input file
juicer_size_file = "hg19" # or "chr6.size"
juicer_hic_file = "GM12878_chr6_25kb_matrix.hic" # output file
juicerBin = " /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar"
stopifnot(file.exists(juicer_pre_file))
arrowhead_TADs_outFile <- "arrowhead_final_domains.txt" 


##### 1) prepare the hic file

juicer_pre_command <- paste("java -Xmx2g -jar",  juicerBin, "pre",
                                  "-n",  # no normalization
                                  "-d",  # only intrachromo
                                  "-r", bin_size,  # resolution
                                  "-c", chromo,  # chromo
                                  juicer_pre_file, juicer_hic_file, # input or output size
                                  juicer_size_file) # size with chr. size or genome (e.g. "hg19")

cat(paste0(juicer_pre_command, "\n"))
system(juicer_pre_command)
# for additional optional arguments: https://github.com/aidenlab/juicer/wiki/Pre

##### 2) run arrowhead
memory_setting <- "-Xms512m -Xmx2048m" # for more memory: "-Xmx20g"

arrowhead_out_folder <- "arrowhead_out"
arrowhead_window <- 2000 # must be even number (Default 2000)

juicer_arrowhead_command <- paste("java", memory_setting, "-jar",  juicerBin, "arrowhead",
                  "-c", chromo, # chromo
                  "-m", arrowhead_window, # sliding window size
                  "-r", bin_size, # resolution
                  "-k NONE ", # normalization
                  juicer_hic_file, arrowhead_out_folder)
cat(paste0(juicer_arrowhead_command, "\n"))
system(juicer_arrowhead_command)

# add --threads for multi-threads

# more info: https://github.com/aidenlab/juicer/wiki/Arrowhead

##### 3) prepare output

arrowhead_out_file <- file.path(arrowhead_out_folder, paste0(bin_size,"_blocks"))
stopifnot(file.exists(arrowhead_out_file))

# x1,x2/y1,y2 = the interval spanned by the domain (contact domains manifest as squares on the diagonal of a Hi-C matrix and as such: x1=y1, x2=y2)
# chr1	x1	x2	chr2	y1	y2	color	score	uVarScore	lVarScore	upSign	loSign
# chr6	8160000	10840000	chr6	8160000	10840000	255,255,0	0.5512261959064583	0.3655323713020236	0.37742302034879593	0.47145328719723184	0.4126297577854671
# chr6	26160000	29600000	chr6	26160000	29600000	255,255,0	0.6551502025255771	0.34722135136955884	0.38974424705160193	0.43868921775898523	0.4260042283298097
# AS I UNDERSTAND THE LAST POSITION GIVES THE LAST POSITION OF THE TAD (NOT THE START OF THE LAST BIN) 

arrowheadData <- read.delim(arrowhead_out_file, header=TRUE, stringsAsFactors = FALSE)
domainsDT_tmp <- arrowheadData[,1:3]  
stopifnot(ncol(domainsDT_tmp) == 3)
colnames(domainsDT_tmp) <- c("chromo", "start", "end")  
stopifnot(is.numeric(domainsDT_tmp$start))
stopifnot(is.numeric(domainsDT_tmp$end))
domainsDT_tmp$start <- domainsDT_tmp$start + 1

# ensure domains are ordered
domainsDT_tmp <- domainsDT_tmp[order(domainsDT_tmp$start,domainsDT_tmp$end ),]
rownames(domainsDT_tmp) <- NULL
domainsDT_tmp <- unique(domainsDT_tmp) # use unique to ensure that I can start >= for start and <= end (both will never be =)
# for each domain -> check if one is nested and if it overlaps
domainsDT <- data.frame(chromo = character(0), start = numeric(0), end = numeric(0), stringsAsFactors = F)

# discard potentially nested or overlapping domains
for(i in seq_len(nrow(domainsDT_tmp))) {
  curr_chromo <- domainsDT_tmp$chromo[i] 
  curr_start <- domainsDT_tmp$start[i] 
  curr_end <- domainsDT_tmp$end[i]
  # can I find a nested domain ? i.e. start between curr_start and curr_end and end between curr_start and curr_end
  if(any(
    domainsDT_tmp$start[-i] >= curr_start & domainsDT_tmp$start[-i] <= curr_end &  
    domainsDT_tmp$end[-i] >= curr_start & domainsDT_tmp$end[-i] <= curr_end 
  )) next
  # if the current domain overlaps a domain already -> just discard it
  if(any(curr_start < domainsDT$end)) next
  stopifnot(curr_start > domainsDT$end & curr_end > domainsDT$end)
  curr_line <- data.frame(chromo = curr_chromo, start = curr_start, end = curr_end, stringsAsFactors = F )
  stopifnot( (curr_end - curr_start + 1) %% binSize == 0)
  domainsDT <- rbind(domainsDT, curr_line)
}
stopifnot(domainsDT$start %in% domainsDT_tmp$start & domainsDT$end %in% domainsDT_tmp$end)
stopifnot(diff(domainsDT$start) > 0)
stopifnot(diff(domainsDT$end) > 0)
stopifnot(domainsDT$end > domainsDT$start)
cat(paste0("... # domains before post-processing\t=\t", nrow(domainsDT_tmp), "\n"))
cat(paste0("... # domains after post-processing\t=\t", nrow(domainsDT), "\n"))

write.table(domainsDT, file=arrowhead_TADs_outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
cat(paste0("... written: ", arrowhead_TADs_outFile))

###########################################################################################################################
# MoC calculation
###########################################################################################################################

source("get_MoC.R")
dt1 <- data.frame(chromo="chr1", start=c(1,101,501,1001), end = c(100,500,1000,2000))
dt2 <- data.frame(chromo="chr1", start=c(1,101,501,1001), end = c(100,500,1000,2000))

get_MoC(dt1,dt2)
get_MoC("moc_test1.txt", "moc_test2.txt")
get_MoC("moc_test1.txt", "moc_test1b.txt")
                    
                    
                    