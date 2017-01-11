#'@include IGSAinput-class.R
setGeneric(name="get_fit", def=function(M, fit_options, use_voom, params) {
    standardGeneric("get_fit")
})

#'@importFrom limma contrasts.fit lmFit treat voom
#'@importFrom stats p.adjust
#'@include ExprData-class.R
#'@include FitOptions-class.R
#'@include SEAparams.R
setMethod(f="get_fit",
    signature=c("ExprData", "FitOptions", "logical", "SEAparams"),
    definition=function(M, fit_options, use_voom, params) {
        act_treat_lfc <- treat_lfc(params);
        act_design <- designMatrix(fit_options);
        act_contrast <- contrast(fit_options);
        act_adj_meth <- adjust_method(params);
        
        if (use_voom) {
            M <- voom(M, design=act_design);
        }
        # Adjust the model
        fit <- lmFit(M, act_design);
        # treat correction
        fit2 <- treat(contrasts.fit(fit, act_contrast), lfc=act_treat_lfc);
        # Adjusted pvalues
        fit2$p.adjust <- apply(fit2$p.value, 2, p.adjust, method=act_adj_meth);
        
        return(fit2);
    }
)

setGeneric(name="igsaGetDEGenes",
    def=function(seaParams, exprData, fitOptions, useVoom) {
    standardGeneric("igsaGetDEGenes")
})

#'@importFrom futile.logger flog.info
#'@include ExprData-class.R
#'@include FitOptions-class.R
#'@include IGSAinput-class.R
setMethod(
    f="igsaGetDEGenes",
    signature=c("SEAparams", "ExprData", "FitOptions", "logical"),
    definition=function(seaParams, exprData, fitOptions, useVoom) {
        act_fit <- get_fit(exprData, fitOptions, useVoom, seaParams);
        
        de_coff <- de_cutoff(seaParams);
        dif <- act_fit$p.adjust[,,drop=FALSE] <= de_coff;
        dif <- unique(rownames(dif)[dif]);
        
        flog.info(paste("DE genes", length(dif), "of a total of", 
                    nrow(exprData), "(", 
                    round(length(dif)/nrow(exprData)*100,2), "%)"));
        
        return(dif);
    }
)
