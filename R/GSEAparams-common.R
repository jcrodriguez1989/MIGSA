#'GSEAparams exploratory functions
#'
#'Several R base overwritten functions to manipulate a GSEAparams object.
#'
#'@param object GSEAparams object.
#'@param ... not in use.
#'
#'@return a summary of the object.
#'
#'@docType methods
#'@name GSEAparams-common
#'@rdname GSEAparams-common
#'
#'@include GSEAparams-class.R
#'@method summary GSEAparams
#'@aliases summary,GSEAparams-method
#'@export summary.GSEAparams
#'@examples
#'gseaParams <- GSEAparams();
#'summary(gseaParams);
#'
summary.GSEAparams <- function(object, ...) {
    stopifnot(validObject(object));
    
    # not showing the other params, as they are from mGSZ, I dont know if 
    # anyone uses them
    res <- object@perm_number;
    names(res) <- "perm_number";
    
    return(res);
}
