# returns the genes present in the IsoDatSet as rows
setMethod(f="dimnames", signature=signature(x="IsoDataSet"),
  definition=function(x) {
    # must return genes that have at least one isoform that is not lowExp
    geneIso <- isoGeneRel(x);
    rown <- as.character(unique(geneIso[,"gene_id"]));
    list(rown, colnames(isoCounts(x)));
})

# number of genes as rows
setMethod(f="dim", signature=signature(x="IsoDataSet"),
  definition=function(x) {
    c(length(rownames(x)), ncol(isoCounts(x)));
})

# returns the IsoDataSet with the isoforms from the genes present in i
setMethod(f="[", signature=c("IsoDataSet", "character"), function(x, i) {
  allGenes <- rownames(x);
  subGenes <- intersect(allGenes, i);
  isoC <- isoCounts(x);
  geneI <- isoGeneRel(x);
  expD <- expData(x);
  
  isoC <- isoC[geneI$gene_id %in% subGenes,];
  geneI <- geneI[geneI$gene_id %in% subGenes,];
  newIsoDataSet <- IsoDataSet(isoC, expD, colnames(expD)[[1]], geneI);
  newIsoDataSet <- buildLowExpIdx(newIsoDataSet);
  return(newIsoDataSet);
})

# returns the IsoDataSet with the isoforms from the genes present in i
setMethod(f="[", signature=c("IsoDataSet", "NULL"), function(x, i) {
  return(NULL);
})
