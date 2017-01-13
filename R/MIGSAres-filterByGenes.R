#'Filters MIGSAres by selected genes
#'
#'\code{filterByGenes} returns a MIGSAres object with only the gene sets which 
#'resulted enriched by at least one gene from a provided list.
#'
#'@param migsaRes MIGSAres object.
#'@param genes character vector of the interest genes for MIGSAres filtering.
#'
#'@return MIGSAres object containing only the gene sets in which genes 
#'contributed to enrichment.
#'
#'@docType methods
#'@name filterByGenes
#'@rdname MIGSAres-filterByGenes
#'
#'@exportMethod filterByGenes
#'
setGeneric(name="filterByGenes", def=function(migsaRes, genes) {
    standardGeneric("filterByGenes")
})

#'@inheritParams filterByGenes
#'@rdname MIGSAres-filterByGenes
#'@aliases filterByGenes,MIGSAres,character-method
#'@include MIGSAres-setEnrCutoff.R
#'@include MIGSAres-genesInSets.R
#'@seealso \code{\link{setEnrCutoff}} and \code{\link{genesInSets}}
#'@examples
#'data(migsaRes);
#'
#'## First lets set a cutoff of 0.1
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.1);
#'
#'## Suppose we are interested in studying these genes:
#'intGenes <- c("g10", "g91", "g388", "g742", "g874");
#'migsaResIntGenes <- filterByGenes(migsaResWCoff, intGenes);
#'
#'## Now in migsaResIntGenes we have the MIGSA results of the gene sets in 
#'## which at least one gene of our list contributed to enrich.
#'migsaResIntGenes
#'
setMethod(
    f="filterByGenes",
    signature=c("MIGSAres", "character"),
    definition=function(migsaRes, genes) {
        stopifnot(validObject(migsaRes));
        stopifnot(length(genes) > 0);
        
        # we must have a cutoff for this function
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        
        # lets get just gene sets enriched in at least one dataset
        actRes <- migsaRes[ rowSums(migsaRes[,-(1:3)], na.rm=!FALSE) > 0, ];
        
        if (!is(actRes, "MIGSAres")) {
            warning("No enriched gene set with used cutOff.");
            return(data.frame());
        }
        migsaRes_genes <- genesInSets(actRes);
        
        if (any(genes %in% colnames(migsaRes_genes))) {
            # if we have any of our interest genes in the enriching ones,
            # then remove the other genes
            migsaRes_genes <- migsaRes_genes[, genes, drop=FALSE];
            
            # keet just the gene sets enriched by at least one of these genes
            migsaRes_genes <- migsaRes_genes[
                        rowSums(migsaRes_genes, na.rm=!FALSE)>0,, drop=FALSE];
            
            # complete gene sets names
            gsNames <- paste(migsaRes[,3,drop=TRUE],
                                migsaRes[,1,drop=TRUE], sep="_");
            
            # so return the subsetted MIGSAres object
            res <- migsaRes[ gsNames %in% rownames(migsaRes_genes),,
                                drop=FALSE];
            
            if (nrow(res) != nrow(migsaRes_genes)) {
                warning("Different number of rows after filtering by genes.");
            }
        } else {
            warning("Provided genes did not contribute to enrich any gene set."
                );
            res <- data.frame();
        }
        
        return(res);
    }
)
