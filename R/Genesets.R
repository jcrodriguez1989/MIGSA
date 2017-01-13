#'@include Genesets-class.R
setGeneric(name="asList", def=function(object) {
    standardGeneric("asList")
})

setMethod(f="asList",
    signature=c("Genesets"),
    definition=function(object) {
        gsets <- object@gene_sets;
        
        # translate a Genesets object to a list. Why? mGSZ uses lists.
        to <- lapply(gsets, function(y) return(genes(y)));
        names(to) <- lapply(gsets, function(y) return(id(y)));
        
        return(to);
    }
)
