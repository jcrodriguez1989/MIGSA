#'MIGSA execution
#'
#'\code{MIGSA} runs a MIGSA execution. Functional analysis is done for each 
#'experiment by means of dEnricher and mGSZ.
#'
#'@param migsaInput MIGSAinput object to analyze.
#'@param ... not in use.
#'@param saveResults logical indicating if back up of each individual 
#'experiment analysis should be saved. It will be saved in 
#'getwd()"/migsaResults/"experimentName".RData".
#'
#'@return MIGSAres object.
#'
#'@docType methods
#'@name MIGSA
#'@rdname MIGSA
#'
#'@exportMethod MIGSA
#'
setGeneric(name="MIGSA", def=function(migsaInput, ...) {
    standardGeneric("MIGSA")
})

#'@inheritParams MIGSA
#'@rdname MIGSA
#'@aliases MIGSA,MIGSAinput,logical-method
#'
#'@importFrom futile.logger flog.info
#'@include IGSA.R
#'@include IGSAinput.R
#'@include IGSAinput-class.R
#'@include IGSAinput-getterSetters.R
#'@include IGSAres.R
#'@include MIGSAres-class.R
#'@include MIGSAinput-getterSetters.R
#'@examples
#'## Lets simulate two expression matrices of 1000 genes and 30 subjects.
#'nGenes <- 1000; # 1000 genes
#'nSamples <- 30; # 30 subjects
#'geneNames <- paste("g", 1:nGenes, sep = ""); # with names g1 ... g1000
#'## Create random gene expression data matrix.
#'set.seed(8818);
#'exprData1 <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
#'rownames(exprData1) <- geneNames;
#'exprData2 <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
#'rownames(exprData2) <- geneNames;
#'
#'## There will be 40 differentialy expressed genes.
#'nDeGenes <- nGenes/25;
#'## Lets generate the offsets to sum to the differentialy expressed genes.
#'deOffsets <- matrix(2*abs(rnorm(nDeGenes*nSamples/2)), ncol=nSamples/2);
#'
#'## Randomly select which are the DE genes.
#'deIndexes1 <- sample(1:nGenes, nDeGenes, replace=FALSE);
#'exprData1[deIndexes1, 1:(nSamples/2)] <-
#'exprData1[deIndexes1, 1:(nSamples/2)] + deOffsets;
#'
#'deIndexes2 <- sample(1:nGenes, nDeGenes, replace=FALSE);
#'exprData2[deIndexes2, 1:(nSamples/2)] <-
#'exprData2[deIndexes2, 1:(nSamples/2)] + deOffsets;
#'
#'exprData1 <- new("MAList",list(M=exprData1));
#'exprData2 <- new("MAList",list(M=exprData2));
#'
#'## 15 subjects with condition C1 and 15 with C2.
#'conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
#'
#'nGSets <- 200; # 200 gene sets
#'## Lets create randomly 200 gene sets, of 10 genes each
#'gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
#'names(gSets) <- paste("set", as.character(1:nGSets), sep="");
#'myGSs <- as.Genesets(gSets, name="myGeneSets");
#'
#'fitOpts <- FitOptions(conditions);
#'
#'igsaInput1 <- IGSAinput(name="igsaInput1", expr_data=exprData1, 
#'fit_options=fitOpts, gene_sets_list=list(myGSs));
#'igsaInput2 <- IGSAinput(name="igsaInput2", expr_data=exprData2, 
#'fit_options=fitOpts, gene_sets_list=list(myGSs));
#'
#'migsaInput <- MIGSAinput(experiments=list(igsaInput1, igsaInput2));
#'
#'## Finally run MIGSA!
#'migsaRes <- MIGSA(migsaInput);
#'
setMethod(
    f="MIGSA",
    signature=c("MIGSAinput"),
    definition=function(migsaInput, saveResults=FALSE) {
        flog.info("*************************************");
        flog.info("Starting MIGSA analysis.");
        
        bp_param <- bpParam(migsaInput);
        geneSets <- geneSetsList(migsaInput);
        
        # for each IGSAinput
        actRes <- lapply(experiments(migsaInput), function(igsaInput) {
            # if migsaInput had gene sets then these must be used
            if (length(geneSets) > 0) {
                geneSetsList(igsaInput) <- geneSets;
            }
            
            igsaRes <- try({ IGSA(igsaInput, bp_param); });
            if (inherits(igsaRes, 'try-error')) return(NA);
            
            # if intermediate results must be saved
            if (saveResults) {
                dir.create("migsaResults", showWarnings=FALSE);
                save(igsaRes, file=paste("migsaResults/", getName(igsaRes),
                    ".RData", sep=""));
            }
            return(igsaRes);
        });
        
        # delete results which gave errors
        actRes <- actRes[!is.na(unlist(actRes))];
        
        migsaRes <- NA;
        if (length(actRes) > 0) {
            # if we have any result then create the MIGSAres
            migsaRes <- MIGSAres(actRes);
        }
        
        return(migsaRes);
    }
)
