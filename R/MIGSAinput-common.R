#'MIGSAinput exploratory functions
#'
#'Several R base overwritten functions to manipulate a MIGSAinput object.
#'
#'@param object MIGSAinput object.
#'@param ... not in use.
#'
#'@return a summary of the object.
#'
#'@docType methods
#'@name MIGSAinput-common
#'@rdname MIGSAinput-common
#'
#'@include IGSAinput.R
#'@include MIGSAinput-class.R
#'@method summary MIGSAinput
#'@aliases summary,MIGSAinput-method
#'@export summary.MIGSAinput
#'@examples
#'## Lets create a basic MIGSAinput object with one IGSAinput object.
#'
#'## First create an expression matrix.
#'maData <- matrix(rnorm(10000),ncol=4);
#'rownames(maData) <- 1:nrow(maData); # It must have rownames (gene names).
#'maExprData <- new("MAList",list(M=maData));
#'
#'## Now lets create the FitOptions object.
#'myFOpts <- FitOptions(c("Cond1", "Cond1", "Cond2", "Cond2"));
#'
#'## Lets create our IGSAinput ready for MIGSA.
#'igsaInput <- IGSAinput(name="igsaInput", expr_data=maExprData, 
#'fit_options=myFOpts);
#'
#'## Finally lets create out MIGSAinput object.
#'migsaInput <- MIGSAinput(experiments=list(igsaInput));
#'summary(migsaInput);
#'
summary.MIGSAinput <- function(object, ...) {
    validObject(object);
    
    # rbind the summary of each IGSAinput
    res <- do.call(rbind, lapply(object@experiments,
                                function(igsaInput) summary(igsaInput)));
    return(res);
}
