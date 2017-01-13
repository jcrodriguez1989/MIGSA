#'Get the genes that contributed to enrichment
#'
#'\code{genesInSets} returns a data.frame with gene sets as rows, genes as 
#'columns, and as value the number of experiments in which each gene 
#'contributed to enrich each gene set. 
#'If it was enriched only by SEA then it returns the genes that contributed in 
#'SEA. If it was enriched only by GSEA then it returns the genes that 
#'contributed in GSEA. If it was enriched by both then it returns the genes 
#'that contributed in both.
#'
#'@param migsaRes MIGSAres object.
#'
#'@return data.frame with the number of experiments in which each gene 
#'contributed to enrich each gene set.
#'
#'@docType methods
#'@name genesInSets
#'@rdname MIGSAres-genesInSets
#'
#'@exportMethod genesInSets
setGeneric(name="genesInSets", def=function(migsaRes) {
    standardGeneric("genesInSets")
})

#'@inheritParams genesInSets
#'@rdname MIGSAres-genesInSets
#'@aliases genesInSets,MIGSAres
#'@include MIGSAres-setEnrCutoff.R
#'@seealso \code{\link{setEnrCutoff}}
#'@examples
#'data(migsaRes);
#'
#'## First lets set a cutoff of 0.1
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.1);
#'gInSets <- genesInSets(migsaResWCoff);
#'class(gInSets); # matrix
#'
#'## Now we can do stuff as check which genes enriched a gene set in all (two) 
#'## experiments.
#'gInSets[rowSums(gInSets==2) > 0, colSums(gInSets==2) > 0];
#'
setMethod(
    f="genesInSets",
    signature=c("MIGSAres"),
    definition=function(migsaRes) {
        stopifnot(validObject(migsaRes));
        
        data_frame_all <- migsaRes@migsa_res_all;
        setkey(data_frame_all, gene_set_name, id);
        
        # we must have a cutoff for this function
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        enr_cutoff <- enrCutoff(migsaRes);
        
        # when cutoff is 1 then it must enrich everything.
        enr_cutoff <- ifelse(enr_cutoff == 1, 1.1, enr_cutoff);
        
        # get the gene sets name and id
        gsets <- apply(unique(data_frame_all[,list(gene_set_name, id)]),2,
                        as.character);
        res <- lapply(split(gsets, seq(nrow(gsets))), function(gset) {
            # for each gene set get the data
            actData <- data_frame_all[ gene_set_name == gset[[1]] &
                                        id == gset[[2]], ]
            enr_genes <- apply(actData, 1, function(y) {
                # for each gene set data
                gsea_enr_genes <- NULL;
                sea_enr_genes <- NULL;
                
                gsea_pval <- as.numeric(y[["GSEA_pval"]]);
                sea_pval  <- as.numeric(y[["SEA_pval"]]);
                
                if (!is.na(gsea_pval) && gsea_pval < enr_cutoff) {
                    # if was enriched by gsea get its enriching genes
                    gsea_enr_genes <- gsub(" ", "", unlist(
                        strsplit(as.character(y[["GSEA_enriching_genes"]]),
                                ",")));
                }
                
                if (!is.na(sea_pval) && sea_pval < enr_cutoff) {
                    # if was enriched by sea get its enriching genes
                    sea_enr_genes <- gsub(" ", "", unlist(
                        strsplit(as.character(y[["SEA_enriching_genes"]]),
                                ",")));
                }
                
                enr_genes <- union(gsea_enr_genes, sea_enr_genes);
                return(enr_genes);
            })
            
            if (is.list(enr_genes)) {
                enr_genes <- do.call(c, enr_genes);
            } else {
                enr_genes <- c(enr_genes);
            }
            
            res <- NA;
            if (length(enr_genes) > 0) {
                # if there was any gene then sort them
                res <- sort(table(enr_genes), decreasing=!FALSE);
            }
            
            return(res);
        })
        names(res) <- paste(gsets[,1], gsets[,2], sep="_");
        
        gsets <- unique(unlist(lapply(res, names)));
        # lets turn these results (list) into a data frame
        finalRes <- do.call(rbind, lapply(names(res), function(y) {
            actRes <- res[[y]];
            newRes <- rep(0, length(gsets));
            if (!is.na(actRes[[1]])) {
                names(newRes) <- gsets;
                newRes[ names(actRes) ] <- actRes;
            }
            return(newRes);
        }))
        rownames(finalRes) <- names(res);
        
        return(finalRes);
    }
)
