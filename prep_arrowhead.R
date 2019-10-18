 # need the size file (white spaces separated ?)
  infile_size <- paste0(outFolder, "/arrowhead_chrsize.txt")
  sizeDT <- data.frame(chr=chromo, size = chrSize)
  cat(paste0("......... write ", infile_size, "\n"))
  # MUST BE TAB SEPARATED !!!
  write.table(sizeDT, file = infile_size, quote=F, row.names=F, col.names=F, sep="\t")
  for(normM in currNorm) {
    # then name for hic file
    outFile <- sub("(.*)\\.(.*?)$", paste0("\\1", "_", normM , ".hic"), inFile )
    outfile_hic <- paste0(outFolder, "/arrowhead_", basename(outFile))
    outFile <- sub("(.*)\\.(.*?)$", paste0("\\1", "_", normM, "_pre.txt"), inFile )
    infile_pre <- paste0(outFolder, "/arrowhead_", basename(outFile))
    
    # whitespace separated file
    # FAQ of Pre juicer -> if already binned/normalized, used "Short with score format"
    # <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <score>
    # pre throws away reads that map to the same restriction fragment. If you use dummy numbers for the frag field, 
    # be sure they are different for the different read ends; that is, <frag1> should be 0 and <frag2> should be 1.
    # readname and strand are also not currently stored within .hic files
    # take the 
    # rao data reindex /binSize + 1 -> so it is zero-based?
    # so I think that if I have interaction chr1 1 40000 & chr1 40001 80000 -> pos1 should be 0 and pos2 should be 40000
    # take similar to sparse matrix 
    if(normM == "LGF") {
      countM <- as.matrix(intdata(htc_object_LGF))
    } else if(normM == "ICE") {
      countM <- as.matrix(intdata(htc_object_ICE))
    }
    sparseMatrix <- Matrix(as.matrix(countM), sparse=T)
    matrixDF <- data.frame(summary(sparseMatrix))
    stopifnot(ncol(matrixDF) == 3)
    # the strand can be both 0 everywhere, but not the fragment
    nRows <- nrow(matrixDF)
    fragVec1 <- rep(0, nRows)
    fragVec2 <- rep(1, nRows)
    posVec1 <- (matrixDF[,1] - 1 ) * binSize
    posVec2 <- (matrixDF[,2] - 1 ) * binSize
    countVec <- matrixDF[,3]
    # chrVec <- rep(sub("chr", "", chromo), nRows)
    chrVec <- rep(chromo, nRows)
    preDT <- data.frame(str1=fragVec1, chr1=chrVec, pos1=posVec1, frag1=fragVec1, 
                        str1=fragVec1, chr1=chrVec, pos1=posVec2, frag1=fragVec2, counts=countVec)
    cat(paste0("......... write ", infile_pre, "\n"))
    write.table(preDT, file=infile_pre, row.names = F, col.names=F, sep=" ", quote=F)
    cat(paste0("......... run juicer tools pre \n"))
    # add -n to not normalize
    # -r binSize -> to do only for a single resolution
    # -c chromo -> only for the chromo
    # -d -> only intrachromo
    command <- paste("java -Xmx2g -jar", path_to_juicer_tools, "pre -n -d -r", binSize, "-c", chromo, infile_pre, outfile_hic, infile_size)
    cat(paste0(command, "\n"))
    system(command)
  }
}
        arrowhead_command="java -Xms512m -Xmx2048m -jar $arrowhead_bin -c $chromo -m 2000 -r $binSize -k NONE $matFile $outFold"
#YUANLONG VERSION FOR MORE MEMORY: java -Xmx20g
#java -Xms512m -Xmx2048m -jar $juicerBin arrowhead -c $chromo -m 2000 -r $binSize -k NONE $hicFile $arrowheadFolder
