#'Lists the gene sets present at enrichr database.
#'
#'\code{enrichrGeneSets} lists the database names present at enrichr.
#'
#'@param pattern character indicating a pattern to filter the database names.
#'
#'@return character with present database names.
#'
#'@docType methods
#'@name enrichrGeneSets
#'@rdname Genesets-enrichrGeneSets
#'
#'@exportMethod enrichrGeneSets
setGeneric(name="enrichrGeneSets", def=function(pattern=".*") {
    standardGeneric("enrichrGeneSets")
})

#'@rdname Genesets-enrichrGeneSets
#'@inheritParams enrichrGeneSets
#'@aliases enrichrGeneSets,ANY-method
#'
#'@importFrom Biobase testBioCConnection
#'@importFrom RJSONIO fromJSON
#'@examples
#'## Lets list all the gene sets that can be downloaded from Enichr website.
#'enrichrGeneSets();
#'
#'## Now lets list only the gene sets that have KEGG in their names.
#'enrichrGeneSets("KEGG");
#'
setMethod(
    f="enrichrGeneSets",
    signature=character(),
    definition=function(pattern=".*") {
        if (!testBioCConnection()) stop("You must have internet connection.");
        
        # enrichr url
        datasetStatisticsUrl <-
                    "http://amp.pharm.mssm.edu/Enrichr/datasetStatistics";
        
        # donwload genesets list (basic info)
        datasetStatistics <- do.call(rbind,
                        fromJSON(datasetStatisticsUrl)$statistics);
        datasetStatistics <- apply(datasetStatistics, 2, unlist);
        datasetStatistics <- data.frame(datasetStatistics);
        
        # giving some format
        datasetStatistics[,2:4] <- apply(datasetStatistics[,2:4], 2,
                                                                as.numeric);
        datasetStatistics[,c(1,5)] <- apply(datasetStatistics[,c(1,5)], 2,
                                                                as.character);
        
        # filtering with pattern (by default returns all)
        datasetStatistics <- datasetStatistics[
            grep(pattern, datasetStatistics[,1], ignore.case=!FALSE), ];
        
        # order by gene set name
        datasetStatistics <- datasetStatistics[order(datasetStatistics[,1]),];
        
        return(datasetStatistics);
    }
)
