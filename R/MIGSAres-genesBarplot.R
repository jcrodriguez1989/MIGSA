# to avoid R CMD check errors we set them as NULL
number = NULL;

#'Genes barplot of the number of gene sets in which each enrich
#'
#'\code{genesBarplot} generates a barplot of the number of gene sets in which 
#'each gene contributed to enrich. Each gene set counts 1 regardless if it was 
#'enriched in many experiments. x-axis each gene, y-axis number of gene sets 
#'in which it contributed to enrich.
#'
#'@param migsaRes MIGSAres object.
#'@param enrFilter numeric. Keep gene sets enriched in at least enrFilter 
#'experiments.
#'@param gsFilter numeric. Keep genes enriched in at least gsFilter gene sets.
#'@param ... not in use.
#'
#'@return ggplot object used as graphic.
#'
#'@docType methods
#'@name genesBarplot
#'@rdname MIGSAres-genesBarplot
#'
#'@exportMethod genesBarplot
setGeneric(name="genesBarplot", def=function(migsaRes, ...) {
    standardGeneric("genesBarplot")
})

#'@inheritParams genesBarplot
#'@rdname MIGSAres-genesBarplot
#'@aliases genesBarplot,MIGSAres-method
#'
#'@importFrom ggplot2 aes geom_bar ggplot
#'@include MIGSAres-genesInSets.R
#'@include MIGSAres-setEnrCutoff.R
#'@seealso \code{\link{setEnrCutoff}}
#'@examples
#'data(migsaRes);
#'
#'## First lets set a cutoff of 0.01
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.01);
#'
#'## Lets check what genes contributed to enrich the highest number of gene 
#'## sets (in more than one gene set).
#'genesBarplot(migsaResWCoff, gsFilter=1);
#'
#'## Moreover we can keep gene sets which where enriched in more than 
#'## enrFilter experiments. To do this, we can use the enrFilter parameter.
#'
setMethod(
    f="genesBarplot",
    signature=c("MIGSAres"),
    definition=function(migsaRes, enrFilter=0, gsFilter=0) {
        stopifnot(validObject(migsaRes));
        
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        
        # keep gene sets enriched in more than enrFilter experiments
        actRes <- migsaRes[ rowSums(migsaRes[,-(1:3)], na.rm=!FALSE) >
                                    enrFilter, ];
        plotGenes <- genesInSets(actRes);
        
        # keep genes enriched in more than gsFilter gene sets
        plotGenes <- plotGenes[, colSums(plotGenes > 0) > gsFilter];
        
        # count 1 disregard if the gene enriched in more than 1 experiment
        plotGenes <- data.frame(id=colnames(plotGenes),
                                number=colSums(plotGenes > 0));
        ## todo: add symbol to plotGenes
        # plotGenes$symbol <- entrez2symbol(plotGenes$id);
        
        # reordering bars
        plotGenes$id <- factor(plotGenes$id,
                                levels=plotGenes$id[order(plotGenes$number,
                                                        decreasing=!FALSE)]);
        plotGenes <- plotGenes[order(plotGenes$number, decreasing=!FALSE),];
        
        p <- ggplot(plotGenes);
        p <- p + geom_bar(aes(id, y=number), stat="identity");
        print(p);
        
        return(p);
    }
)
