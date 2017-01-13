#'Gets the Gene Ontology depths of specified ids
#'
#'\code{getHeights} returns the heights of given a list of ids (GO IDs).
#'
#'@param ids character vector indicating the queried GO ids.
#'@param minHeight logical indicating if the minimum or maximum height must be 
#'calculated. If it is FALSE then the longest path to the root is calculated, 
#'otherwise, the shortest path.
#'@param ... not in use.
#'
#'@return list with each term height.
#'
#'@docType methods
#'@name getHeights
#'@rdname MIGSAres-getHeights
#'
#'@exportMethod getHeights
#'@examples
#'## Lets load breast cancer results.
#'data(bcMigsaRes);
#'
#'## And get the first 40 Gene Ontology gene sets IDs.
#'goIds <- bcMigsaRes[bcMigsaRes$GS_Name %in% c("BP", "CC", "MF"), "id"];
#'fst40goIds <- goIds[1:40,];
#'
#'## And lets get the heights in the GO tree structure.
#'getHeights(fst40goIds);
#'
setGeneric(name="getHeights", def=function(ids, ...) {
    standardGeneric("getHeights")
})

#'@inheritParams getHeights
#'@rdname MIGSAres-getHeights
#'@aliases getHeights,character-method
#'
#'@importFrom GO.db GOBPPARENTS GOCCPARENTS GOMFPARENTS
#'@include GoAnalysis.R
#'
setMethod(
    f="getHeights",
    signature=c("character"),
    definition=function(ids, minHeight=TRUE) {
        allParents <- as.list(GOBPPARENTS);
        allParents <- c(allParents, as.list(GOMFPARENTS));
        allParents <- c(allParents, as.list(GOCCPARENTS));
        
        # for each GO id find its height
        result <- lapply(ids, function(x) {
            getHeight(x, minHeight, allParents)
        });
        
        return(unlist(result));
    }
)
