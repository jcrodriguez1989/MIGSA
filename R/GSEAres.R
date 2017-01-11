#'@include GenesetRes.R
GSEAres <- setClass(
    Class="GSEAres",
    slots=c(
        gene_sets_res="list"
    ),
    prototype=list(
    ),
    validity=function(object) {
        gene_sets_res_ok <- all(unlist(lapply(object@gene_sets_res, function(x)
            is(x, "GenesetRes"))));

        return(gene_sets_res_ok);
    }
)

#'@importFrom data.table as.data.table data.table setkey
setMethod(
    f="as.data.table",
    signature=c("GSEAres"),
    definition=function(x, wGenesInfo=FALSE, ...) {
        to <- do.call(rbind, lapply(x@gene_sets_res, function(res) {
            actInfo <- asCharacter(res, wGenesInfo=wGenesInfo);
            return(actInfo);
        }))
        to <- data.table(to);
        if (wGenesInfo) {
            colnames(to) <- c("id", "name", "GSEA_enriched", "GSEA_score",
                                "GSEA_pval", "GSEA_enriching_genes",
                                "GSEA_GS_genes");
        } else {
            colnames(to) <- c("id", "name", "GSEA_enriched", "GSEA_score",
                                "GSEA_pval");
        }
        setkey(to, id, name);

        return(to)
    }
)

setGeneric(name="gene_sets_res", def=function(object) {
    standardGeneric("gene_sets_res")
})

setMethod(f="gene_sets_res", signature="GSEAres",
    definition=function(object) {
        return(object@gene_sets_res)
    }
)
