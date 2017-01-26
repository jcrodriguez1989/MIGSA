#'Creates a GeneSetCollection object using the Gene Ontology data base
#'
#'\code{loadGo} creates a GeneSetCollection object from the data present at 
#'the Gene Ontology data base (org.Hs.eg.db R package).
#'
#'@param ontology character indicating which ontology must be loaded. 
#'Must be one of BP, MF or CC.
#'
#'@return A GeneSetCollection object  (Genes are with their EntrezGene ID).
#'
#'@docType methods
#'@name loadGo
#'@rdname Genesets-loadGo
#'@seealso \code{\link{as.Genesets}}
#'@seealso \code{\link{Genesets-enrichr}}
#'@seealso \code{\link{geneSetsFromFile}}
#'
#'@exportMethod loadGo
setGeneric(name="loadGo", def=function(ontology) {
    standardGeneric("loadGo")
})

#'@rdname Genesets-loadGo
#'@inheritParams loadGo
#'@aliases loadGo,character-method
#'
#'@importFrom AnnotationDbi as.list Ontology Term
#'@importFrom GSEABase GOCollection GeneSet EntrezIdentifier GeneSetCollection
#'@importFrom org.Hs.eg.db org.Hs.egGO2ALLEGS
#'@importFrom futile.logger flog.error
#'@examples
#'## Lets load the Cellular Components gene sets from the Gene Ontology.
#'\dontrun{
#'ccGSets <- loadGo("CC");
#'}
# biocCheck gives me a note if I dont have any line ran.
#'\dontshow{
#'aux <- 1:10; rm(aux);
#'}
#'
setMethod(
    f="loadGo",
    signature=c("character"),
    definition=function(ontology) {
        ontology <- toupper(ontology);
        if (!ontology %in% c("BP", "CC", "MF")) {
            stop("Ontology must be one of: 'BP', 'CC', 'MF'.");
        }
        
        # download GO gene sets
        go <- org.Hs.egGO2ALLEGS;
        go <- as.list(go);
        go <- lapply(go, unique);
        
        goCollection <- GOCollection(ontology=ontology,
            evidenceCode=as.character(NA));
        
        # filter the desired ontology
        go <- go[Ontology(names(go)) == ontology];
        gsets <- lapply(names(go), function(x) {
            gset <- GeneSet(go[[x]], setIdentifier=Term(x), setName=x,
                geneIdType=EntrezIdentifier(), 
                collectionType=goCollection);
            
            return(gset);
        });
        
        gsetsColl <- GeneSetCollection(gsets);
        return(gsetsColl);
    }
)
