#'Creates a Genesets from file
#'
#'\code{geneSetsFromFile} creates a Genesets object from the data present in 
#'a file. The file must be a tab separated values file (tsv). Each line will 
#'parse to a Geneset. First field will be the Geneset id, the second the name 
#'and the remaining are the genes.
#'
#'@param filePath character with the path of the file to parse.
#'@param name character with the name of this Genesets.
#'@param is_GO logical indicating if this Genesets are from the Gene Ontology. 
#'If true, then each Geneset id must be a GO id.
#'@param ... not in use.
#'
#'@return Genesets object.
#'
#'@docType methods
#'@name geneSetsFromFile
#'@rdname Genesets-geneSetsFromFile
#'
#'@include Genesets-class.R
#'@include Geneset.R
#'@exportMethod geneSetsFromFile
#'@examples
#'## Create some fake gene sets in a data.frame to save them in disk and then
#'## load them (10 gene sets with 20 genes each (it is not neccesary that they
#'## have the same number of genes).
#'gsets <- data.frame(
#'     IDs=paste("set", 1:10), 
#'     Names=rep("", 10), 
#'     matrix(paste("gene", 1:(10*20)), nrow=10));
#'
#'## And save this file as a tab separated file.
#'write.table(gsets, file="fakeGsets.tsv", sep="\t",
#'     col.names=FALSE, row.names=FALSE, quote=FALSE);
#'
#'## Now lets load this tsv file as a Genesets object.
#'myGsets <- geneSetsFromFile("fakeGsets.tsv", "fakeGsets");
#'
#'## And lets delete this tsv file (so we dont have garbage in our disk).
#'unlink("fakeGsets.tsv");
#'
setGeneric(name="geneSetsFromFile", def=function(filePath, name, ... ) {
    standardGeneric("geneSetsFromFile")
})

#'@rdname Genesets-geneSetsFromFile
#'@inheritParams geneSetsFromFile
#'@aliases geneSetsFromFile,character-method
#'
#'@importFrom futile.logger flog.error
setMethod(
    f="geneSetsFromFile",
    signature=c("character"),
    definition=function(filePath, name, is_GO=FALSE) {
        stopifnot(length(filePath) == 1);
        
        # check if file exists, its not a directory, its readable
        if (!file.exists(filePath) || dir.exists(filePath)
                            || (file.access(filePath, 4) == -1)) {
            flog.error("Gene sets file path error, check if it is readable.");
        }
        
        tmp <- readLines(filePath);
        gsets <- lapply(tmp, function(actLine) {
            # for each line, fst item is id, snd is name, rest are genes
            t <- strsplit(actLine,'\t')[[1]]
            actGS <- Geneset(id=t[[1]], name=t[[2]],
                                genes=unique(t[3:length(t)]));
        })
        
        return(Genesets(name=name, is_GO=is_GO, gene_sets=gsets));
    }
)
