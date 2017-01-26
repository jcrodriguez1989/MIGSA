#'IGSAinput S4 class implementation in R
#' 
#'This S4 class contains all the necessary inputs to execute a functional 
#'analysis (SEA and GSEA) on one experiment.
#'Important: Make sure that gene IDs are concordant between the expression 
#'matrix and the provided gene sets.
#'
#'@slot name character indicating the name of this experiment.
#'@slot expr_data ExprData object with the expression data (MicroArray or 
#'RNAseq).
#'@slot use_voom logical indicating wether to use 
#'\code{\link[limma]{voom}} before fit or not. If using RNAseq data it is 
#'recommended to set it to TRUE (default: FALSE).
#'@slot fit_options FitOptions object with the parameters to be used when 
#'fitting the model.
#'@slot gene_sets_list named list of GeneSetCollection objects to be 
#'tested for enrichment (names must be unique).
#'@slot sea_params SEAparams object with the parameters to be used 
#'by SEA.
#'@slot gsea_params GSEAparams object with the parameters to be used 
#'by GSEA.
#'
#'@docType methods
#'@name IGSAinput-class
#'@rdname IGSAinput-class
#'@seealso \code{\link{ExprData-class}}
#'@seealso \code{\link{SEAparams-class}}
#'@seealso \code{\link{GSEAparams-class}}
#'@seealso \code{\link{IGSAinput-getterSetters}}
#'@seealso \code{\link{getDEGenes}}
#'@seealso \code{\link{MIGSA}}
#'@seealso \code{\link{summary}}
#'
#'@importFrom futile.logger flog.error
#'@importFrom GSEABase GeneSet GeneSetCollection
#'@include ExprData-class.R
#'@include FitOptions-class.R
#'@include SEAparams.R
#'@include GSEAparams.R
#'@include Genesets.R
#'@export IGSAinput
#'@examples
#'## Lets create a basic IGSAinput object.
#'## First create a expression matrix.
#'maData <- matrix(rnorm(10000),ncol=4);
#'rownames(maData) <- 1:nrow(maData); # It must have rownames (gene names).
#'maExprData <- new("MAList",list(M=maData));
#'
#'## Now lets create the FitOptions object.
#'myFOpts <- FitOptions(c("Cond1", "Cond1", "Cond2", "Cond2"));
#'
#'## Finally lets create the Genesets to test for enrichment.
#'library(GSEABase);
#'myGs1 <- GeneSet(as.character(1:10), setIdentifier="fakeId1", 
#'     setName="fakeName1");
#'myGs2 <- GeneSet(as.character(7:15), setIdentifier="fakeId2", 
#'     setName="fakeName2");
#'myGSs <- GeneSetCollection(list(myGs1, myGs2));
#'
#'## And now we can create our IGSAinput ready for MIGSA.
#'igsaInput <- IGSAinput(name="igsaInput", expr_data=maExprData, 
#'fit_options=myFOpts, gene_sets_list=list(myGSs=myGSs));
#'
IGSAinput <- setClass(
    Class="IGSAinput",
    slots=c(
        name="character",
        expr_data="ExprData",
        use_voom="logical",
        fit_options="FitOptions",
        gene_sets_list="list",
        sea_params="SEAparams",
        gsea_params="GSEAparams"
    ),
    prototype=list(
        use_voom=FALSE
    ),
    validity=function(object) {
        # must have name
        name_ok <- length(object@name) == 1 && object@name != "";
        
        # must have genes and samples
        expr_data_ok <- (!is.null(object@expr_data)) &&
                            ncol(object@expr_data) > 1 &&
                            nrow(object@expr_data) > 1;
        
        # gene_sets_list is a list of GeneSetCollection
        gene_sets_list_ok <- all(unlist(lapply(object@gene_sets_list,
                                    function(x) is(x, "GeneSetCollection"))));
        
        # GeneSetCollection names must be unique
        if (gene_sets_list_ok) {
            gss_names <- names(object@gene_sets_list);
            gene_sets_list_ok <- length(gss_names) ==
                                        length(unique(gss_names));
            
            # every gene set collection must have a name
            gene_sets_list_ok <- gene_sets_list_ok && 
                length(gss_names) == length(object@gene_sets_list);
            
            if (!gene_sets_list_ok) {
                flog.error("GeneSetCollection names must be unique");
            }
        }
        
        # check that the FitOptions and the ExprData are concordant
        fit_opts_ok <- nrow(designMatrix(object@fit_options)) ==
                                                    ncol(object@expr_data);
        
        return(name_ok && expr_data_ok && gene_sets_list_ok && fit_opts_ok);
    }
)
