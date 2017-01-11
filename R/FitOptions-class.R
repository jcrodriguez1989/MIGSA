#'FitOptions S4 class implementation in R
#' 
#'This S4 class contains the parameters to provide for model fitting.
#'See \link{FitOptions} to create one FitOptions object.
#'
#'@docType methods
#'@name FitOptions-class
#'@rdname FitOptions-class
#'
#'@importFrom futile.logger flog.error
#'@importFrom stats model.matrix
#'@exportClass FitOptions
setClass(
    Class="FitOptions",
    slots=c(
        col_data="data.frame",
        formula="formula",
        contrast="numeric",
        design_matrix="matrix"
    ),
    prototype=list(
    ),
    validity=function(object) {
        design_matrix <- object@design_matrix;
        formula <- object@formula;
        col_data <- object@col_data;
        contrast <- object@contrast;
        
#         design_matrix_error <- all.equal(design_matrix, model.matrix(formula,
#                                             data=col_data));
#         contst_names <- make.names(
#                     colnames(model.matrix(formula, data=col_data)));
        contrast_ok <- length(contrast) == ncol(design_matrix);
        if (!contrast_ok) {
            flog.error("Contrast length must be equal to design_matrix
                                                        columns number.");
        }
        
        return(contrast_ok);
    }
)
