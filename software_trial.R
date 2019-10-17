###########################################################################################################################
# TopDom - package version
###########################################################################################################################


topDom_file <- "/mnt/etemp/marie/TADcall_yuanlong/25kb/input_caller/chr6/GM12878_chr6_25kb_matrix_pos_zero.txt"

topDom_outFile <- "topDom_results"
topDom_TADs_outFile <- "topDom_final_domains.txt"

if(!require(TopDom)) {
  if( !require(remotes)) {
    install.packages("remotes")
    library(remotes)
  }
  remotes::install_github("HenrikBengtsson/TopDom")
  library(topDom)
}

  
topDom_out <- TopDom(topDom_file, window.size=5, outFile=topDom_outFile)
  # > head(topDom_out[["domain"]])
  # chr from.id from.coord to.id to.coord    tag   size
  # 1  chr6       1          0     6   150000    gap 150000
  # 9  chr6       7     150000    16   400000 domain 250000
  # 10 chr6      17     400000    23   575000 domain 175000
  # 11 chr6      24     575000    32   800000 domain 225000
  # 12 chr6      33     800000    57  1425000 domain 625000
  # 13 chr6      58    1425000    64  1600000 domain 175000
  # > head(topDom_out[["bed"]])                                                                                                                                                                                 
  # chrom chromStart chromEnd   name
  # 1   chr6          0   150000    gap
  # 9   chr6     150000   400000 domain
  # 10  chr6     400000   575000 domain
  # 11  chr6     575000   800000 domain
  # 12  chr6     800000  1425000 domain
  # 13  chr6    1425000  1600000 domain
  # > head(topDom_out[["binSignal"]])                                                                                                                                                                           
  # id  chr from.coord to.coord local.ext  mean.cf pvalue
  # 1  1 chr6          0    25000      -0.5   0.0000      1
  # 2  2 chr6      25000    50000      -0.5   0.0000      1
  # 3  3 chr6      50000    75000      -0.5   0.0000      1
  # 4  4 chr6      75000   100000      -0.5 392.5379      1
  # 5  5 chr6     100000   125000      -0.5 271.3408      1
  
topDom_tads_dt <- topDom_out[["bed"]]
topDom_tads_dt <- topDom_tads_dt[as.character(topDom_tads_dt$name) == "domain",]
out_dt <- topDom_tads_dt[,c("chrom", "chromStart","chromEnd")]

write.table(out_dt, file=topDom_TADs_outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
cat(paste0("... written: ", topDom_TADs_outFile))

  

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
# once installed
library(CaTCH)

catch_file <- "/mnt/etemp/marie/TADcall_yuanlong/25kb/input_caller/chr6/GM12878_chr6_25kb_matrix_list.txt"

catch_file <- "/mnt/etemp/marie/TADcall_yuanlong/25kb/input_caller/chr6/GM12878_chr6_25kb_matrix_list.txt"


ri_thresh <- 0.65


# domain.call(input)
# domain.call.parallel(inputs,ncpu=parallel::detectCores()-1L)
# 
# domain.call	
# A list with two elements
# ncluster: A data.frame with 3 columns
# -Chromosome (chromosome)
# -Reciprocal insulation (RI)
# -Number of domains (ndomains)
# clusters: A data.frame with 4 columns
# -Chromosome (chromosome)
# -Reciprocal insulation (RI)
# -Start of domain (start)
# -End of domain (end)
# domain.call.parallel	
# A list of outputs as in domain.call for each file


CaTCH_TADs_outFile