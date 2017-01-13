#'IGSAinput exploratory functions
#'
#'Several R base overwritten functions to manipulate a IGSAinput object.
#'
#'@param object IGSAinput object.
#'@param ... not in use.
#'
#'@return a summary of the object.
#'
#'@docType methods
#'@name IGSAinput-common
#'@rdname IGSAinput-common
#'
#'@include IGSAinput-class.R
#'@method summary IGSAinput
#'@aliases summary,IGSAinput-method
#'@export summary.IGSAinput
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
#'## And now we can create our IGSAinput ready for MIGSA.
#'igsaInput <- IGSAinput(name="myIgsaInput", expr_data=maExprData, 
#'fit_options=myFOpts);
#'summary(igsaInput);
#'
summary.IGSAinput <- function(object, ...) {
    validObject(object);
    
    deGenes <- getDEGenes(object);
    # number of samples of each contrast
    ctrst <- table(col_data(deGenes@fit_options));
    sea_params <- summary(deGenes@sea_params);
    gsea_params <- summary(deGenes@gsea_params);
    
    res <- c(deGenes@name,
            ncol(deGenes@expr_data),
            paste(names(ctrst), collapse="VS"),
            ctrst[[1]],
            ctrst[[2]],
            length(deGenes@gene_sets_list),
            deGenes@use_voom,
            nrow(deGenes@expr_data),
            sea_params,
            gsea_params,
            round(100*(as.numeric(sea_params[[4]]) / 
                nrow(deGenes@expr_data)), 2)
    );
    names(res) <- c("exp_name",
                    "#samples",
                    "contrast",
                    "#C1",
                    "#C2",
                    "#gene_sets",
                    "use_voom",
                    "#genes",
                    names(sea_params),
                    names(gsea_params),
                    "%de_genes");
    return(res);
}
