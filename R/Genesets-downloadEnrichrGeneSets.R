#'Downloads Genesets from enrichr database
#'
#'\code{downloadEnrichrGeneSets} creates a list of Genesets objects 
#'downloading the specified ones from enrichr website 
#'(http://amp.pharm.mssm.edu/Enrichr/).
#'
#'@param geneSetNames list of characters with the names of the gene sets to 
#'download. Must be listed at \code{\link{enrichrGeneSets}}.
#'@param deleteMultipleEntrez logical indicating if multiple Entrez IDs should 
#'be deleted or repeated.
#'Note: Enrichr uses Gene Symbol, if org.Hs.eg translates it into several 
#'Entrez IDs then if deleteMultipleEntrez == FALSE all Entrez are removed else 
#'all are included.
#'@param ... not in use.
#'
#'@return list of Genesets objects.
#'
#'@docType methods
#'@name downloadEnrichrGeneSets
#'@rdname Genesets-downloadEnrichrGeneSets
#'
#'@include Genesets-class.R
#'@include Genesets-enrichrGeneSets.R
#'@include Geneset.R
#'@exportMethod downloadEnrichrGeneSets
setGeneric(name="downloadEnrichrGeneSets", def=function(geneSetNames, ...) {
    standardGeneric("downloadEnrichrGeneSets")
})

#'@rdname Genesets-downloadEnrichrGeneSets
#'@inheritParams downloadEnrichrGeneSets
#'@aliases downloadEnrichrGeneSets,character-method
#'
#'@importFrom AnnotationDbi as.list
#'@importFrom Biobase testBioCConnection
#'@importFrom futile.logger flog.info
#'@importFrom org.Hs.eg.db org.Hs.egALIAS2EG
#'@examples
#'## Lets download BioCarta gene sets from Enrichr.
#'## Make sure you use the same names as listed with enrichrGeneSets() .
#'bioCartaGSs <- downloadEnrichrGeneSets(c("BioCarta_2015"));
#'
setMethod(
    f="downloadEnrichrGeneSets",
    signature=c("character"),
    definition=function(geneSetNames,
    deleteMultipleEntrez=!FALSE) {
        
        downloadUrlFst <-
"http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=";
        libraries <- enrichrGeneSets()[,1];
        
        geneSetNames <- intersect(geneSetNames, libraries);
        if (length(geneSetNames) < 1) {
            stop("No gene sets found to download from Enrichr");
        }
        
        if (!testBioCConnection()) stop("You must have internet connection.");
        
        symbol2entrez <- as.list(org.Hs.egALIAS2EG);
        
        libraries <- lapply(geneSetNames, function(libName) {
            actUrl <- paste(downloadUrlFst, libName, sep="");
            tmp <- readLines(actUrl);
            flog.info(paste("Downloaded", libName));
            tmp <- strsplit(tmp, "\t");
            flog.info("Converting Symbol to Entrez");
            geneSets <- lapply(tmp, function(actGS) {
                actId <- actGS[[1]];
                actName <- actGS[[2]];
                entrez <- lapply(actGS[2:length(actGS)], function(symbol) {
                    translates <- symbol2entrez[[symbol]];
                    if (length(translates) > 1 && deleteMultipleEntrez) {
                        translates <- NULL;
                    }
                    
                    return(translates);
                })
                entrez <- Reduce(union, entrez);
                
                res <- NA;
                
                if (length(entrez) > 0) {
                    res <- Geneset(id=actId, name=actName, genes=entrez);
                }
                
                return(res);
            })
            
            geneSets <- geneSets[!is.na(geneSets)];
            res <- NA;
            is_GO <- grepl("^GO_", libName);
            if (length(geneSets) > 0) {
                res <- Genesets(name=libName, gene_sets=geneSets, is_GO=is_GO);
            }
            
            return(res);
        })
        names(libraries) <- geneSetNames;
        return(libraries);
    }
)
