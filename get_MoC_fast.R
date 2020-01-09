
# yuanlong version
library('GenomicRanges')

  get_MoC_GenomicRanges = function(set1DT, set2DT)
  {
    colnames(set1DT) = colnames(set2DT) = c('chr', 'start', 'end')
    refGR = makeGRangesFromDataFrame(set1DT)
    testGR = makeGRangesFromDataFrame(set2DT)

    hits <- findOverlaps(refGR, testGR)
    overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
    # refGR[queryHits(hits)]
    Q = width(testGR[subjectHits(hits)])
    P = width(refGR[queryHits(hits)])
    F = width(overlaps)

    data = data.frame(P=P, Q=Q, F=F)
    MoC_score = sum(F^2/P/Q)

    MoC_score <- MoC_score/(sqrt(nrow(set1DT) * nrow(set2DT)) - 1) 
    return(MoC_score)
  }

