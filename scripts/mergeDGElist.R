mergeDGEList <- function(firstDgeList, secondDgeList,
                         DGEListLabels=NULL) {
  stopifnot(is.null(DGEListLabels) || length(DGEListLabels)==2)
  commonFeat <- intersect(rownames(firstDgeList$counts),
                          rownames(secondDgeList$counts))
  commonSampleCols <- intersect(colnames(firstDgeList$samples),
                                colnames(secondDgeList$samples))
  
  mergedCounts <- cbind(firstDgeList$counts[commonFeat,],
                        secondDgeList$counts[commonFeat,])
  mergedGenes <- firstDgeList$genes[commonFeat,]
  mergedSamples <- rbind(firstDgeList$samples[, commonSampleCols],
                         secondDgeList$samples[, commonSampleCols])
  if(!is.null(DGEListLabels) & length(DGEListLabels)==2) {
    labels <- factor(rep(DGEListLabels,
                         c(nrow(firstDgeList$samples),
                           nrow(secondDgeList$samples))),
                     levels=DGEListLabels)
    mergedSamples$DGEListLabel <- labels
  }
  
  uniqSampleNames <- make.unique(colnames(mergedCounts))
  colnames(mergedCounts) <- rownames(mergedSamples) <- uniqSampleNames
  
  res <- DGEList(counts=mergedCounts,
                 samples=mergedSamples,
                 genes=mergedGenes)
  # res$samples <- ribiosUtils::removeColumns(res$samples,
  #                                           c("lib.size.1", "norm.factors.1"))
  return(res)
}
