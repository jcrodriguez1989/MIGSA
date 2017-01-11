#'FitOptions constructor
#'
#'\code{FitOptions} creates a FitOptions object. If the vector of samples is
#'provided (must be two different) then it will contrast Condition1 vs. 
#'Condition2. If not, it should be provided with a data.frame x, the formula 
#'and the contrast, it will create the model matrix using x as data, and the 
#'formula.
#'
#'@param x There are two options for x:
#'\itemize{
#'\item It can be a character vector containing the two conditions (length 
#'must be the same as the number of subjects to use).
#'\item It can be a data.frame used as data by 
#'\code{\link[stats]{model.matrix}}.
#'}
#'@param formula (only used if x is data.frame) used by 
#'\code{\link[stats]{model.matrix}}.
#'@param contrast (only used if x is data.frame) the contrast to test.
#'@param ... not in use.
#'
#'@return FitOptions object.
#'
#'@include FitOptions-class.R
#'@docType methods
#'@name FitOptions
#'@rdname FitOptions
#'@export FitOptions
FitOptions <- function(x, ...) {
    UseMethod("FitOptions", x);
}

#'@rdname FitOptions
#'@inheritParams FitOptions
#'@aliases FitOptions.default
#'@importFrom stats model.matrix
#'@export FitOptions.default
#'@examples
#'## Supose we have 15 subjects, the first 8 from Condition1 and the last 7 
#'## from Condition2, lets create the corresponding FitOptions object to test
#'## Condition1 vs. Condition2.
#'l <- c(rep("Condition1", 8), rep("Condition2", 7));
#'fit_options <- FitOptions(l);
#'
FitOptions.default <- function(x, ...) {
    if (length(x) < 2) {
        stop("More than two labels required.");
    }
    if (length(unique(x)) != 2) {
        stop("Exactly two possible conditions required.");
    }
    
    act_col_data <- data.frame(cond=factor(x));
    act_formula  <- ~cond-1;
    
    act_contrast <- c(-1, 1);
    .Object <- FitOptions(x=act_col_data, formula=act_formula,
                            contrast=act_contrast);
    return(.Object);
}

#'@rdname FitOptions
#'@inheritParams FitOptions
#'@aliases FitOptions.data.frame
#'@importFrom stats model.matrix
#'@export FitOptions.data.frame
#'@examples
#'## Otherwise if we have the data and formula for model.matrix function and 
#'## the desired contrast, we can create the FitOptions object as:
#'myData <- data.frame(cond=c(rep("Condition1", 8), rep("Condition2", 7)));
#'myFormula <- ~cond - 1;
#'myContrast <- c(-1, 1);
#'fit_options <- FitOptions(myData, myFormula, myContrast);
#'
FitOptions.data.frame <- function(x, formula, contrast, ...) {
    stopifnot(is(x, "data.frame"));
    stopifnot(is(formula, "formula"));
    stopifnot(is(contrast, "numeric"));

    act_design <- model.matrix(formula, data=x);
    .Object <- new("FitOptions", col_data=x, formula=formula, 
                    contrast=contrast, design_matrix=act_design);
    return(.Object);
}
