#'@include GenesetRes.R
SEAres <- setClass(
    Class="SEAres",
    slots=c(
        gene_sets_res="list"
    ),
    prototype=list(
    ),
    validity=function(object) {
        gene_sets_res_ok <- all(unlist(lapply(object@gene_sets_res,
            function(x) is(x, "GenesetRes"))));
        
        return(gene_sets_res_ok);
    }
)

#'@importFrom data.table as.data.table data.table setkey
as.data.table.SEAres <- function(x, wGenesInfo=FALSE, ...) {
    # convert each gene set result to character vector
    to <- do.call(rbind, lapply(x@gene_sets_res, function(res) {
        actInfo <- asCharacter(res, wGenesInfo=wGenesInfo);
        return(actInfo);
    }))
    
    # make it a data.table and put the colname
    to <- data.table(to);
    if (wGenesInfo) {
        colnames(to) <- c("id", "name", "SEA_enriched", "SEA_score",
                        "SEA_pval", "SEA_enriching_genes", "SEA_GS_genes");
    } else {
        colnames(to) <- c("id", "name", "SEA_enriched", "SEA_score",
                        "SEA_pval");
    }
    # set data.table keys, to merge faster
    setkey(to, id, name);
    
    return(to)
}

#'@include GSEAres.R
setMethod(f="gene_sets_res", signature="SEAres",
    definition=function(object) {
        return(object@gene_sets_res)
    }
)
