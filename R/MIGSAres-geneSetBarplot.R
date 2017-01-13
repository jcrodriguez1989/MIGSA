# to avoid R CMD check errors we set them as NULL
number = NULL;

#'Gene set barplot
#'Gene set barplot of the number of experiments in which each enrich
#'
#'\code{geneSetBarplot} generates a barplot of the number of experiments in 
#'which each gene set was enriched. x-axis each gene set, y-axis times it was 
#'enriched (0 to #experiments).
#'
#'@param migsaRes MIGSAres object.
#'@param enrFilter numeric. Keep gene sets enriched in at least enrFilter 
#'experiments.
#'@param ... not in use.
#'
#'@return ggplot object used as graphic.
#'
#'@docType methods
#'@name geneSetBarplot
#'@rdname MIGSAres-geneSetBarplot
#'
#'@exportMethod geneSetBarplot
setGeneric(name="geneSetBarplot", def=function(migsaRes, ...) {
    standardGeneric("geneSetBarplot")
})

#'@inheritParams geneSetBarplot
#'@rdname MIGSAres-geneSetBarplot
#'@aliases geneSetBarplot,MIGSAres-method
#'
#'@importFrom ggplot2 aes geom_bar ggplot
#'@include MIGSAres-setEnrCutoff.R
#'@seealso \code{\link{setEnrCutoff}}
#'@examples
#'data(migsaRes);
#'
#'## First lets set a cutoff of 0.1
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.1);
#'
#'## Lets check in how many experiments each gene set was enriched (in more 
#'## than one experiment).
#'geneSetBarplot(migsaResWCoff, enrFilter=1);
#'
setMethod(
    f="geneSetBarplot",
    signature=c("MIGSAres"),
    definition=function(migsaRes, enrFilter=0) {
        stopifnot(validObject(migsaRes));
        
        # we must have a cutoff for this function
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        
        # keep gene sets enriched in more than enrFilter experiments
        actRes <- migsaRes[ rowSums(migsaRes[,-(1:3)], na.rm=!FALSE) >
                                enrFilter, ];
        plotRes <- data.frame(actRes[,1:3],
                                number=rowSums(actRes[,-(1:3)], na.rm=!FALSE));
        
        # reordering bars
        plotRes$id <- factor(plotRes$id,
                            levels=plotRes$id[order(plotRes$number,
                                                    decreasing=!FALSE)]);
        plotRes <- plotRes[order(plotRes$number, decreasing=!FALSE),];
        
        p <- ggplot(plotRes);
        p <- p + geom_bar(aes(id, y=number), stat="identity");
        print(p);
        
        return(p);
    }
)
