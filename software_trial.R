options(scipen=100)

startTime <- Sys.time()

# ***************************************************************************************************************************************************
# # *************************** PART 0: input data *************************************************************************************************
# ***************************************************************************************************************************************************


# download data from: 
# https://raw.githubusercontent.com/CSOgroup/TAD-benchmarking-scripts/master/Example_Data/GM12878_chr19_25kb_matrix_pos_zero.txt.gz
# unzip the file, save as:
input_file <- "GM12878_chr19_25kb_matrix_pos_zero.txt"
chromo <- "chr19"
binKb <- 25
bin_size <- binKb*1000

topDom_file <- input_file

# ***************************************************************************************************************************************************
# # *************************** PART I: TAD calling *************************************************************************************************
# ***************************************************************************************************************************************************

runTopDom <- FALSE
runCaTCH <- FALSE
runArrowhead <- FALSE
runHiCseg <- FALSE

runMoC <- TRUE
runPlotMoC <- FALSE

if(runTopDom){
cat(paste0("> START TopDom\n"))
###########################################################################################################################
# TopDom - package version
###########################################################################################################################
### package installation:
# - in R console: require(remotes); remotes::install_github("HenrikBengtsson/TopDom") # might need to install.packages("remotes")

stopifnot(file.exists(topDom_file))
topDom_outFile <- file.path(paste0("out_TopDom_", chromo), "TopDom_results")
topDom_TADs_outFile <- file.path(paste0("out_TopDom_", chromo), "TopDom_final_domains.txt")
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
# system.time(topDom_out <- TopDom(topDom_file, window.size=5, outFile=topDom_outFile))
# user  system elapsed 
# 54.278   0.632  54.915 

##### 2) prepare output
topDom_tads_dt <- topDom_out[["bed"]]
topDom_tads_dt <- topDom_tads_dt[as.character(topDom_tads_dt$name) == "domain",]
domainsDT <- topDom_tads_dt[,c("chrom", "chromStart","chromEnd")]
# to 1-based start positions:
domainsDT$chromStart <- domainsDT$chromStart + 1 
# ensure ordering
domainsDT <- domainsDT[order(domainsDT$chromStart, domainsDT$chromEnd),]
stopifnot(domainsDT$end > domainsDT$start)
stopifnot(diff(domainsDT$end) > 0)
stopifnot(diff(domainsDT$start) > 0)
write.table(domainsDT, file=topDom_TADs_outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
cat(paste0("... ", nrow(domainsDT) , " domains written in: ", topDom_TADs_outFile, "\n"))
}
###########################################################################################################################
# TopDom - original script version
###########################################################################################################################

# original site not longer available
# script still available: https://github.com/HenrikBengtsson/TopDom/blob/0.0.2/R/TopDom.R

source("TopDom.R")
# => then use of TopDom() function similar as above

if(runCaTCH){
cat(paste0("> START CaTCH\n"))
###########################################################################################################################
# CaTCH
###########################################################################################################################

### package installation:
# - download CaTCH_1.0.tar.gz from https://github.com/zhanyinx/CaTCH_R
# - in a terminal: R CMD check CaTCH_1.0.tar.gz # to check that you have all libraries
# - in a terminal: R CMD INSTALL CaTCH_1.0.tar.gz # to install the package

source("infile_convert.R")
catch_file <- paste0("GM12878_", chromo, "_", binKb, "_", "kb_matrix_list.txt")
stopifnot(file.exists(topDom_file))
TopDom_to_CaTCH(infile=topDom_file, outfile=catch_file, outZeroBased=FALSE)
stopifnot(file.exists(catch_file))

# installation instructions: https://github.com/zhanyinx/CaTCH_R/
# once installed
stopifnot(require(CaTCH))
stopifnot(file.exists(catch_file))
catch_TADs_outFile <- file.path(paste0("out_CaTCH_", chromo), "CaTCH_final_domains.txt")
dir.create(dirname(catch_TADs_outFile), recursive = TRUE)

ri_thresh <- 0.65
bin_size <- 25*1000

##### 1) run CaTCH
CaTCH_out <- domain.call(catch_file)
# system.time(CaTCH_out <- domain.call(catch_file))
# user  system elapsed 
# 35.072   1.337  36.501 


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
ref_ri <- all_ri[which.min(abs(all_ri - ri_thresh))]  # not exact, e.g. 0.65 is infact 0.6499999999...
# ref_ri = 0.7 -> 306 domains
# ref_ri = 0.65 -> 389 domains
# ref_ri = 0.6 -> 526 domains
cat(paste0("... selected RI thresh\t=\t", round(ri_thresh,4),"\n"))
cat(paste0("... best matching RI thresh\t=\t", round(ref_ri,4),"\n"))
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
stopifnot(domainsDT$end > domainsDT$start)
stopifnot(diff(domainsDT$end) > 0)
stopifnot(diff(domainsDT$start) > 0)

write.table(domainsDT, file=catch_TADs_outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
cat(paste0("... ", nrow(domainsDT) , " domains written in: ", catch_TADs_outFile), "\n")


# multi-threads (if > 1 chromo) domain.call.parallel(inputs,ncpu=parallel::detectCores()-1L)
}

if(runArrowhead){
cat(paste0("> START Arrowhead\n"))
###########################################################################################################################
# arrowhead
###########################################################################################################################
# installation instruction: https://github.com/aidenlab/juicer/wiki/Installation

source("infile_convert.R")
  
juicer_pre_file <- paste0("GM12878_", chromo, "_25kb_matrix.pre") # input file
stopifnot(file.exists(topDom_file))
TopDom_to_arrowhead(infile=topDom_file,
                    outfile=juicer_pre_file)
stopifnot(file.exists(juicer_pre_file))


juicer_size_file <- paste0(chromo, ".size") # hg19 # ! chr6.size should be tab-separated ! [for using hg19 might be 6, not chr6]
juicer_hic_file <- paste0("GM12878_", chromo, "_", binKb, "kb_matrix.hic") # output file
juicerBin <- " /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar"
arrowhead_TADs_outFile <- file.path(paste0("out_Arrowhead_", chromo), "arrowhead_final_domains.txt")
dir.create(dirname(arrowhead_TADs_outFile), recursive = TRUE)


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
# system.time(system(juicer_pre_command))
# user  system elapsed 
# 376.311   5.578 379.448 
# for additional optional arguments: https://github.com/aidenlab/juicer/wiki/Pre

##### 2) run arrowhead
memory_setting <- "-Xms512m -Xmx2048m" # for more memory: "-Xmx20g"
arrowhead_out_folder <- paste0("out_Arrowhead_", chromo)
arrowhead_window <- 2000 # must be even number (Default 2000)

juicer_arrowhead_command <- paste("java", memory_setting, "-jar",  juicerBin, "arrowhead",
                  "-c", chromo, # chromo
                  "-m", arrowhead_window, # sliding window size
                  "-r", bin_size, # resolution
                  "-k NONE ", # normalization
                  juicer_hic_file, arrowhead_out_folder)
cat(paste0(juicer_arrowhead_command, "\n"))
system(juicer_arrowhead_command)
# system.time(system(juicer_arrowhead_command))
# Arrowhead complete
# user  system elapsed 
# 33.669   3.213  34.543 


# add --threads for multi-threads
# more info: https://github.com/aidenlab/juicer/wiki/Arrowhead

##### 3) prepare output
# (retain smallest level of the hierarchy and non-overlapping)
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

# discard if potentially nested or overlapping domains
for(i in seq_len(nrow(domainsDT_tmp))) {
  curr_chromo <- domainsDT_tmp$chromo[i] 
  curr_start <- domainsDT_tmp$start[i] 
  curr_end <- domainsDT_tmp$end[i]
  # can I find a nested domain ? i.e. start between curr_start and curr_end and end between curr_start and curr_end
  # if yes -> discard (I retain smallest level of hierarchy)
  if(any(
    domainsDT_tmp$start[-i] >= curr_start & domainsDT_tmp$start[-i] <= curr_end &  
    domainsDT_tmp$end[-i] >= curr_start & domainsDT_tmp$end[-i] <= curr_end 
  )) next
  # if the current domain overlaps a domain already -> just discard it
  if(any(curr_start < domainsDT$end)) next
  stopifnot(curr_start > domainsDT$end & curr_end > domainsDT$end)
  curr_line <- data.frame(chromo = curr_chromo, start = curr_start, end = curr_end, stringsAsFactors = F )
  stopifnot( (curr_end - curr_start + 1) %% bin_size == 0)
  domainsDT <- rbind(domainsDT, curr_line)
}
stopifnot(domainsDT$start %in% domainsDT_tmp$start & domainsDT$end %in% domainsDT_tmp$end)
stopifnot(diff(domainsDT$start) > 0)
stopifnot(diff(domainsDT$end) > 0)
stopifnot(domainsDT$end > domainsDT$start)
cat(paste0("... # domains before post-processing\t=\t", nrow(domainsDT_tmp), "\n"))
cat(paste0("... # domains after post-processing\t=\t", nrow(domainsDT), "\n"))

# discdDomains <- which(! domainsDT_tmp$start %in% domainsDT$start)
# domainsDT_tmp[1:3,]
# ensure ordering
domainsDT <- domainsDT[order(domainsDT$start, domainsDT$end),]
stopifnot(domainsDT$end > domainsDT$start)
stopifnot(diff(domainsDT$end) > 0)
stopifnot(diff(domainsDT$start) > 0)

write.table(domainsDT, file=arrowhead_TADs_outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
cat(paste0("... ", nrow(domainsDT) , " domains written in: ", arrowhead_TADs_outFile, "\n"))


}
if(runHiCseg){
cat(paste0("> START HiCseg\n"))
###########################################################################################################################
# HiCseg
###########################################################################################################################

hicseg_TADs_outFile <- file.path(paste0("out_HiCseg_", chromo), "HiCseg_final_domains.txt")
dir.create(dirname(hicseg_TADs_outFile), recursive = TRUE)

# ! # of change points
hicseg_model <- "D"
hicseg_dist <- "G" # Gaussian for normalized data
hicseg_maxChangePoints <- 1000 # might be slow if too many

if(!require(HiCseg)) {
  install.packages("HiCseg")
  stopifnot(require(HiCseg))
}


hicseg_dt <- read.delim(topDom_file, header=FALSE)
stopifnot(ncol(hicseg_dt) == nrow(hicseg_dt) + 3)
hicseg_dt <- hicseg_dt[,-c(1:3)]
stopifnot(ncol(hicseg_dt) == nrow(hicseg_dt))
hicseg_mat <- as.matrix(hicseg_dt)
stopifnot(ncol(hicseg_mat) == nrow(hicseg_mat))

mat_size <- dim(hicseg_mat)[1] #6843

##### 1) run HiCseg (! may be slow !)

HiCseg_out <- HiCseg_linkC_R(mat_data = hicseg_mat,
                            size_mat = mat_size, # Size of the data matrix
                             nb_change_max = hicseg_maxChangePoints, # Maximal number of change-points
                             distrib = hicseg_dist, # "B" is for Negative Binomial distribution, "P" is for the Poisson distribution and "G" is for the Gaussian distribution.
                             model = hicseg_model) #  "D" for block-diagonal and "Dplus" for the extended block-diagonal model.

# system.time():
#     user   system  elapsed 
# 1861.582    6.838 1868.378 
  

# t_hat	: Contains the estimated change-points
# J	: Values of the log-likelihood for different number of change-points up to some constants
# t_est_mat	 :It gives the matrix of the estimated change-points for different number of change-points
# head(HiCseg_out[["t_hat"]])
# [1]  3  5  7  9 13 16
# head(HiCseg_out[["J"]])
# [1] -569707711663 -566712112603 -565114308231 -563626723790 -562195217366 -560876084365
# head(HiCseg_out[["t_est_mat"]])

### 2) prepare output

# => from the vignette, if the 1st value is 39, draws the domain border at 40 -> the 39th position is included in the 1st domain
# => hence, if the first valuees are 2 and 4, for 25kb data, the domains should range 1-50'000 and 50'001-100'000

changePoints_bin <- HiCseg_out[["t_hat"]]
changePoints_bin <- changePoints_bin[changePoints_bin > 0]  # because if it found less change points than nb_change_max, wil fill with 0

chromoSize <- mat_size * bin_size
changePoints <- changePoints_bin * bin_size
curr_start <- 1
domainsDT <- data.frame(chromo=character(0), start=numeric(0), end=numeric(0))  

for(point in changePoints) {
  stopifnot( ((point - curr_start + 1) %% bin_size) == 0 )
  tmpDT <- data.frame(chromo = chromo, start=curr_start, end = point)
  domainsDT <- rbind(domainsDT, tmpDT)
  curr_start <- point + 1
}
if(curr_start < chromoSize){
  lastDT <- data.frame(chromo = chromo, start = curr_start, end = chromoSize)
  domainsDT <- rbind(domainsDT, lastDT)
}
# ensure ordering
domainsDT <- domainsDT[order(domainsDT$start, domainsDT$end),]

stopifnot(domainsDT$end > domainsDT$start)
stopifnot(diff(domainsDT$end) > 0)
stopifnot(diff(domainsDT$start) > 0)

write.table(domainsDT, file=hicseg_TADs_outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
cat(paste0("... ", nrow(domainsDT) , " domains written in: ", hicseg_TADs_outFile, "\n"))
}





# ***************************************************************************************************************************************************
# # *************************** PART II: comparison of genomic partitions ***************************************************************************
# ***************************************************************************************************************************************************



if(runMoC){
cat(paste0("> START calculate MoC \n"))
###########################################################################################################################
# MoC calculation
###########################################################################################################################
source("get_MoC.R")
  
topdom_dt <- read.delim(file.path(paste0("out_TopDom_", chromo), "TopDom_final_domains.txt"), header=FALSE, col.names=c("chromo", "start", "end"))
catch_dt <- read.delim(file.path(paste0("out_CaTCH_", chromo), "CaTCH_final_domains.txt"), header=FALSE, col.names=c("chromo", "start", "end"))
arrowhead_dt <- read.delim(file.path(paste0("out_Arrowhead_", chromo), "arrowhead_final_domains.txt"), header=FALSE, col.names=c("chromo", "start", "end"))
hicseg_dt <- read.delim(file.path(paste0("out_HiCseg_", chromo), "HiCseg_final_domains.txt"), header=FALSE, col.names=c("chromo", "start", "end"))

cat("> MoC TopDom vs. Arrowhead:\n")
cat(get_MoC(topdom_dt, arrowhead_dt), "\n") # 0.54
cat("> MoC TopDom vs. CaTCH:\n")
cat(get_MoC(topdom_dt, catch_dt), "\n") # 0.65
cat("> MoC  TopDom vs. HiCseg:\n")
cat(get_MoC(topdom_dt, hicseg_dt), "\n") # 0.69

cat("> MoC Arrowhead vs. CaTCH:\n")
cat(get_MoC(arrowhead_dt, catch_dt), "\n") # 0.55
cat("> MoC Arrowhead vs. HiCseg:\n")
cat(get_MoC(arrowhead_dt, hicseg_dt), "\n") # 0.42

cat("> MoC CaTCH vs. HiCseg:\n")
cat(get_MoC(catch_dt, hicseg_dt), "\n") # 0.57

}
  
if(runPlotMoC){
cat(paste0("> START plot MoC heatmap\n"))
    
###########################################################################################################################
# MoC calculation and similarity heatmap
###########################################################################################################################

source("plot_TADlist_similarities.R")
  
outHeightGG <- 7
outWidthGG <- 7
plotFile <- paste0(chromo, "_moc_plot.svg")

topdom_dt <- read.delim(file.path(paste0("out_TopDom_", chromo), "TopDom_final_domains.txt"), header=FALSE, col.names=c("chromo", "start", "end"))
catch_dt <- read.delim(file.path(paste0("out_CaTCH_", chromo), "CaTCH_final_domains.txt"), header=FALSE, col.names=c("chromo", "start", "end"))
arrowhead_dt <- read.delim(file.path(paste0("out_Arrowhead_", chromo), "arrowhead_final_domains.txt"), header=FALSE, col.names=c("chromo", "start", "end"))
hicseg_dt <- read.delim(file.path(paste0("out_HiCseg_", chromo), "HiCseg_final_domains.txt"), header=FALSE, col.names=c("chromo", "start", "end"))



moc_plot <- plot_TADlist_comparison(TAD_list= list(TopDom=topdom_dt, CaTCH=catch_dt, Arrowhead=arrowhead_dt, HiCseg=hicseg_dt),
                                    plotTit = paste0("Comparison of TAD calling with MoC"),
                                    "get_MoC", nCpu=1)

ggsave(filename = plotFile, plot = moc_plot[[1]], height=outHeightGG, width=outWidthGG)
cat(paste0("... MoC heatmap written in: ", plotFile, "\n"))


}

# 
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################

cat("*** DONE\n")
cat(paste0(startTime, " - ", Sys.time(), "\n"))

