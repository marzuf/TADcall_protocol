chromo <- "chr19"
TopDom_dt <- read.delim(file.path(paste0("out_TopDom_", chromo), "TopDom_final_domains.txt"),
                          col.names = c("chromo", "start", "end"),
                          header=F, stringsAsFactors = FALSE)
TopDom_dt$region <- paste0("chr19_TAD", 1:nrow(TopDom_dt), "_TopDom")

CaTCH_dt <- read.delim(file.path(paste0("out_CaTCH_", chromo), "CaTCH_final_domains.txt"),
col.names = c("chromo", "start", "end"),
header=F, stringsAsFactors = FALSE)
CaTCH_dt$region <- paste0("chr19_TAD", 1:nrow(CaTCH_dt), "_CaTCH")

Arrowhead_dt <- read.delim(file.path(paste0("out_Arrowhead_", chromo), "arrowhead_final_domains.txt"),
col.names = c("chromo", "start", "end"),
header=F, stringsAsFactors = FALSE)
Arrowhead_dt$region <- paste0("chr19_TAD", 1:nrow(Arrowhead_dt), "_Arrowhead")


HiCseg_dt <- read.delim(file.path(paste0("out_HiCseg_", chromo), "HiCseg_final_domains.txt"),
col.names = c("chromo", "start", "end"),
header=F, stringsAsFactors = FALSE)
HiCseg_dt$region <- paste0("chr19_TAD", 1:nrow(HiCseg_dt), "_HiCseg")

# 2140001
toplot_start <- TopDom_dt$start[10] - 10000
# 2960000
toplot_end <- TopDom_dt$end[13] + 10000

toplot_td <- TopDom_dt[TopDom_dt$start >= toplot_start & TopDom_dt$end <= toplot_end,]
toplot_ar <- Arrowhead_dt[Arrowhead_dt$start >= toplot_start & Arrowhead_dt$end <= toplot_end,]
toplot_hs <- HiCseg_dt[HiCseg_dt$start >= toplot_start & HiCseg_dt$end <= toplot_end,]
toplot_ca <- CaTCH_dt[CaTCH_dt$start >= toplot_start & CaTCH_dt$end <= toplot_end,]

all_to_plot <- do.call(rbind, list(toplot_td, toplot_ar, toplot_hs, toplot_ca))

write.table(all_to_plot, col.names = FALSE, row.names = FALSE, sep="\t", quote=F, file = "select_TADs_to_plot.txt")