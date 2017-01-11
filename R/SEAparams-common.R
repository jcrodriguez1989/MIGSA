#'SEAparams exploratory functions
#'
#'Several R base overwritten functions to manipulate a SEAparams object.
#'
#'@param object SEAparams object.
#'@param ... not in use.
#'
#'@return a summary of the object.
#'
#'@docType methods
#'@name SEAparams-common
#'@rdname SEAparams-common
#'
#'@include SEAparams-class.R
#'@method summary SEAparams
#'@aliases summary,SEAparams-method
#'@export summary.SEAparams
#'@examples
#'seaParams <- SEAparams();
#'summary(seaParams);
#'
summary.SEAparams <- function(object, ...) {
    stopifnot(validObject(object));
    
    br <- object@br;
    if (length(br) > 1) {
        br <- "UserDefined";
    }
    
    res <- c(object@treat_lfc,
            object@de_cutoff,
            object@adjust_method,
            length(object@de_genes),
            br);
    names(res) <- c("treat_lfc", "de_cutoff", "adjust_method", "#de_genes",
                                                                    "br");
    
    return(res);
}
