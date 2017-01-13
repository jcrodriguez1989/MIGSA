#'Genesets S4 class implementation in R
#' 
#'This S4 class represents a collection of gene sets.
#'
#'@slot name character containing the 
#'name of this collection of gene sets.
#'@slot gene_sets list of Geneset, containing each gene set (can not have 
#'repeated Geneset ids).
#'@slot is_GO logical indicating wether this sets are from Gene 
#'Ontology (for further analyses and plots). If true, then each Geneset id 
#'must be a GO id.
#'
#'@importFrom futile.logger flog.error
#'@include Geneset.R
#'@docType methods
#'@name Genesets-class
#'@rdname Genesets-class
#'@export Genesets
#'@examples
#'## Lets create manually a gene sets containing three gene set.
#'## First lets create manually three gene set.
#'myGs1 <- Geneset(id="fakeId1", name="fakeName1", genes=as.character(1:10));
#'myGs2 <- Geneset(id="fakeId2", name="fakeName2", genes=as.character(7:15));
#'myGs3 <- Geneset(id="fakeId3", name="fakeName2", genes=as.character(20:28));
#'
#'## Now we can create our Genesets object.
#'myGSs <- Genesets(name="myGenesets", gene_sets=list(myGs1, myGs2, myGs3),
#'is_GO=FALSE);
#'
Genesets <- setClass(
    Class="Genesets",
    slots=c(
        name="character",
        gene_sets="list",
        is_GO="logical"
    ),
    prototype=list(
        is_GO=FALSE
    ),
    validity=function(object) {
        # name cant be empty
        name_ok <- object@name != "" && length(object@name) == 1;
        
        # gene_sets must be a list of Geneset
        gene_sets_ok <- all(unlist(lapply(object@gene_sets, function(x)
                                                        is(x, "Geneset"))));
        
        # must have at least one Geneset
        gene_sets_ok <- gene_sets_ok && length(object@gene_sets) > 0;
        
        if (gene_sets_ok) {
            # Geneset names must be unique
            gs_ids <- unlist(lapply(object@gene_sets, function(x)
                                                            return(id(x))));
            
            gene_sets_ok <- length(gs_ids) == length(unique(gs_ids));
            if (!gene_sets_ok) {
                flog.error("Into each Genesets, Geneset ids must be unique");
            }
        }
        
        return(name_ok && gene_sets_ok);
    }
)
