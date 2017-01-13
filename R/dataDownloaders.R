#'Functions to download different RDatas
#'
#'Different functions to download datasets used in MIGSA's vignette.
#'
#'@param dataName character name of the data to download, one of
#'\itemize{
#'\item "tcgaMAdata" (~ 8.8 MB) tcga microarray expresion matrix and subtypes.
#'\item "tcgaRNAseqData" (~ 7.8 MB) tcga RNAseq expresion matrix and subtypes.
#'\item "pbcmcData" (~ 51 MB) pbcmc microarray expresion matrices and subtypes.
#'\item "bcMigsaRes" (~ 15 MB) breast cancer MIGSA results.
#'}
#'
#'@return list with downloaded data.
#'
#'@docType methods
#'@name loadMIGSAdata
#'@rdname dataDownloaders
#'
#'@exportMethod loadMIGSAdata
#'@examples
#'## Lets load tcga MicroArray data.
#'tcgaMAdata <- loadMIGSAdata("tcgaMAdata");
#'
#'## Lets load tcga RNAseq data.
#'tcgaRNAseqData <- loadMIGSAdata("tcgaRNAseqData");
#'
#'## Lets load pbcmc microarray datasets.
#'pbcmcData <- loadMIGSAdata("pbcmcData");
#'
#'## Lets load breast cancer MIGSA results.
#'bcMigsaRes <- loadMIGSAdata("bcMigsaRes");
#'
setGeneric(name="loadMIGSAdata", def=function(dataName) {
    standardGeneric("loadMIGSAdata")
})

#'@inheritParams loadMIGSAdata
#'@rdname dataDownloaders
#'@aliases loadMIGSAdata
#'
#'@importFrom Biobase testBioCConnection
#'@importFrom utils download.file
#'
setMethod(
    f="loadMIGSAdata",
    signature=c("character"),
    definition=function(dataName) {
#         if (inDevelopment) return(TRUE);
        
        # checking that just one name is provided and it is valid
        if (length(dataName) > 1 ||
            !dataName %in% c("tcgaMAdata", "tcgaRNAseqData", "pbcmcData", 
                            "bcMigsaRes", "mGSZspeedup", "testing")) {
            stop(paste("dataName must be one of",
                "tcgaMAdata, tcgaRNAseqData, pbcmcData, bcMigsaRes"));
        }
        
        rdata <- paste(dataName, ".RData", sep="");
        
        # I should check if this works in Windows
        tmpFile <- paste(tempdir(), "..", rdata, sep="/");
        
        if (!testBioCConnection()) stop("You must have internet connection.");
        
        # if file was already downloaded then skip this
        if (!file.exists(tmpFile)) {
            # just for testing purposes
            if (dataName == "testing") return(TRUE)
            
            # download the file in temporary path (/tmp/ in unix)
            download.file(
                paste("https://github.com/jcrodriguez1989/MIGSAdata/",
                    "raw/master/", rdata, sep=""), tmpFile);
        }
        res <- get(load(tmpFile));
        
        return(res);
    }
)
