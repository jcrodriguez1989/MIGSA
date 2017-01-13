#'Plots genes heatmap
#'
#'\code{genesHeatmap} plots a heatmap with the number of experiments in which 
#'each gene contributed to enrich each gene set.
#'
#'@param migsaRes MIGSAres object.
#'@param enrFilter numeric. Keep gene sets enriched in at least enrFilter 
#'experiments.
#'@param gsFilter numeric. Keep genes enriched in at least gsFilter gene sets.
#'@param ... additional parameters passed to heatmap.2 function.
#'
#'@return list returned by heatmap.2 function (plotted data).
#'
#'@docType methods
#'@name genesHeatmap
#'@rdname MIGSAres-genesHeatmap
#'
#'@exportMethod genesHeatmap
setGeneric(name="genesHeatmap", def=function(migsaRes, ...) {
    standardGeneric("genesHeatmap")
})

#'@inheritParams genesHeatmap
#'@rdname MIGSAres-genesHeatmap
#'@aliases genesHeatmap,MIGSAres-method
#'
#'@importFrom gplots heatmap.2
#'@include MIGSAres-genesInSets.R
#'@include MIGSAres-setEnrCutoff.R
#'@seealso \code{\link{setEnrCutoff}}
#'@examples
#'data(migsaRes);
#'
#'## First lets set a cutoff of 0.1
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.1);
#'
#'## Lets check what genes contributed to enrich the highest number of gene 
#'## sets (in more than one gene set).
#'genesHeatmap(migsaResWCoff, gsFilter=1);
#'
#'## Moreover we can keep gene sets which where enriched in more than 
#'## enrFilter experiments. To do this, we can use the enrFilter parameter.
#'
setMethod(f="genesHeatmap",
    signature=c("MIGSAres"),
    definition=function(migsaRes, enrFilter=0, gsFilter=0, ...) {
        stopifnot(validObject(migsaRes));
        
        # we must have a cutoff for this function
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        
        # keep gene sets enriched in more than enrFilter experiments
        actRes <- migsaRes[ rowSums(migsaRes[,-(1:3)], na.rm=!FALSE) >
                                    enrFilter, ];
        plotGenes <- genesInSets(actRes);
        
        # keep genes enriched in more than gsFilter gene sets
        plotGenes <- plotGenes[, colSums(plotGenes > 0) > gsFilter];
        
        p <- heatmap.2(plotGenes, labRow=rep("", nrow(plotGenes)), ...);
        p$data <- plotGenes;
        
        return(p);
    }
)
