#'MIGSAinput S4 class implementation in R
#' 
#'This S4 class contains all the necessary inputs to execute a MIGSA.
#'
#'@slot experiments list of IGSAinput objects to execute.
#'@slot gene_sets_list list of Genesets to be tested for enrichment. 
#'If set in MIGSAinput class then it will be tested in every IGSAinput.
#'@slot bp_param (optional) BiocParallelParam to execute MIGSA 
#'in parallel.
#'
#'@docType methods
#'@name MIGSAinput-class
#'@rdname MIGSAinput-class
#'
#'@importClassesFrom BiocParallel BiocParallelParam
#'@importFrom BiocParallel bpparam
#'@importFrom futile.logger flog.error
#'@include IGSAinput.R
#'@export MIGSAinput
#'@examples
#'## Lets create a basic MIGSAinput object with two IGSAinput objects.
#'## These two IGSAinput objects will be identical except for their ExprData.
#'
#'## First create an expression matrix.
#'maData1 <- matrix(rnorm(10000),ncol=4);
#'rownames(maData1) <- 1:nrow(maData1); # It must have rownames (gene names).
#'maExprData1 <- new("MAList",list(M=maData1));
#'
#'maData2 <- matrix(rnorm(10000),ncol=4);
#'rownames(maData2) <- 1:nrow(maData2); # It must have rownames (gene names).
#'maExprData2 <- new("MAList",list(M=maData2));
#'
#'## Now lets create the FitOptions object.
#'myFOpts <- FitOptions(c("Cond1", "Cond1", "Cond2", "Cond2"));
#'
#'## Lets create the Genesets to test for enrichment.
#'myGs1 <- Geneset(id="fakeId1", name="fakeName1", genes=as.character(1:10));
#'myGs2 <- Geneset(id="fakeId2", name="fakeName2", genes=as.character(7:15));
#'myGSs <- Genesets(name="myGenesets", gene_sets=list(myGs1, myGs2),
#'is_GO=FALSE);
#'
#'## Lets create our IGSAinputs ready for MIGSA.
#'myIgsaInput1 <- IGSAinput(name="myIgsaInput1", expr_data=maExprData1, 
#'fit_options=myFOpts, gene_sets_list=list(myGSs));
#'
#'myIgsaInput2 <- IGSAinput(name="myIgsaInput2", expr_data=maExprData2, 
#'fit_options=myFOpts, gene_sets_list=list(myGSs));
#'
#'## Finally lets create out MIGSAinput object.
#'myMigsaInput <- MIGSAinput(experiments=list(myIgsaInput1, myIgsaInput2));
#'
MIGSAinput <- setClass(
    Class="MIGSAinput",
    slots=c(
        experiments="list",
        gene_sets_list="list",
        bp_param="BiocParallelParam"
    ), prototype=list(
        bp_param=bpparam()
    ),
    validity=function(object) {
        # experiments is a list of IGSAinputs
        experiments_ok <- all(unlist(lapply(object@experiments, function(x)
                                                        is(x, "IGSAinput"))));
        experiments_ok <- experiments_ok && length(object@experiments) > 0;
        
        # IGSAinput names must be unique
        if (experiments_ok) {
            exprs_names <- unlist(lapply(object@experiments, function(x)
                                                            return(name(x))));
            experiments_ok <- all(length(exprs_names) ==
                                    length(unique(exprs_names)));
            if (!experiments_ok) {
                flog.error("IGSAinput names must be unique");
            }
        }
        
        # gene_sets_list is a list of Genesets
        gene_sets_list_ok <- all(unlist(lapply(object@gene_sets_list,
                                    function(x) is(x, "Genesets"))));
        
        # Genesets names must be unique
        if (gene_sets_list_ok) {
            gss_names <- unlist(lapply(object@gene_sets_list, function(x)
                                return(name(x))));
            gene_sets_list_ok <- all(length(gss_names) ==
                                        length(unique(gss_names)));
            if (!gene_sets_list_ok) {
                flog.error("Genesets names must be unique");
            }
        }
        
        return(experiments_ok && gene_sets_list_ok);
    }
)
