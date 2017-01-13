#'Creates a Genesets from a list
#'
#'\code{as.Genesets} creates a Genesets object from the data present in a list.
#'Each element will parse to a Geneset. For each list element, its name will be 
#'the Geneset id, and the content are the genes.
#'
#'@param x list of character vectors which are the genes corresponding to each 
#'Geneset. The list must have names (unique).
#'@param name character indicating the name of these gene sets e.g. "KEGG" 
#'(default: "no_name").
#'@param ... not in use.
#'
#'@return Genesets object.
#'
#'@docType methods
#'@name as.Genesets
#'@rdname Genesets-as.Genesets
#'
#'@include Genesets-class.R
#'@include Geneset.R
#'@exportMethod as.Genesets
setGeneric("as.Genesets", def=function(x, ...) {
    standardGeneric("as.Genesets")
})

#'@rdname Genesets-as.Genesets
#'@inheritParams as.Genesets
#'@aliases as.Genesets,list-method
#'@examples
#'## Lets create a list with three manually created gene sets and load it as a 
#'## Genesets object.
#'myGs1 <- as.character(1:10);
#'myGs2 <- as.character(15:21);
#'myGs3 <- as.character(25:30);
#'myGssList <- list(myGs1, myGs2, myGs3);
#'names(myGssList) <- c("myGs1", "myGs2", "myGs3");
#'myGss <- as.Genesets(myGssList, name="myGeneSets");
#'
setMethod(
    f="as.Genesets",
    signature=c("list"),
    definition=function(x, name="no_name") {
        if (length(unique(names(x))) != length(x)) {
            stop("The list must have all unique names.");
        }
        act_gss <- lapply(names(x), function(act_id) {
            # these are the genes
            act_genes <- unique(x[[act_id]]);
            
            # if there is any gene then create the Geneset
            if (length(act_genes) == 0 ||
                (length(act_genes) == 1 && act_genes == "")) {
                act_gs <- NULL;
            } else {
                act_gs <- Geneset(id=act_id, genes=act_genes);
            }
            return(act_gs);
        })
        # delete invalid genesets
        act_gss <- act_gss[!unlist(lapply(act_gss, is.null))];
        
        if (length(act_gss) < 1) {
            stop("Could not create any Geneset object,
                    maybe they where empty");
        }
        res <- Genesets(name=name, gene_sets=act_gss, is_GO=FALSE);
        return(res)
    }
)
