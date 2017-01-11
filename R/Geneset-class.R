#'Geneset S4 class implementation in R
#' 
#'This S4 class represents an individual gene set.
#'
#'@slot id character containing the id of this gene set.
#'@slot name character containing the name of this gene set (can be 
#'empty).
#'@slot genes character vector of the genes conforming this gene set 
#'(must contain at least one gene).
#'
#'@export Geneset
#'@docType methods
#'@name Geneset-class
#'@rdname Geneset-class
#'@examples
#'## Lets create manually a gene set
#'fakeGenes <- c("8818", "2406", "64176");
#'myGs <- Geneset(id="fakeId", name="fakeName", genes=fakeGenes);
#'
Geneset <- setClass(
    Class="Geneset",
    slots=c(
        id="character",
        name="character",
        genes="character"
    ),
    prototype=list(
    ),
    validity=function(object) {
        id_ok <- object@id != "" && length(object@id) == 1; # id cant be empty
        genes_ok <- length(object@genes) > 0; # there is at least one gene
        
        if (genes_ok && length(object@genes) == 1) {
            genes_ok <- object@genes != ""; # there is at least one gene
        }
        
        return(id_ok && genes_ok);
    }
)
