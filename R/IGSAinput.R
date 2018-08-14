#'@include IGSAinput-class.R
setGeneric(name="get_fit", def=function(M, fit_options, params) {
    standardGeneric("get_fit")
})

#'@importFrom limma contrasts.fit lmFit treat voom
#'@importFrom stats p.adjust
#'@include ExprData-class.R
#'@include FitOptions-class.R
#'@include SEAparams.R
setMethod(f="get_fit",
    signature=c("ExprData", "FitOptions", "SEAparams"),
    definition=function(M, fit_options, params) {
        act_treat_lfc <- treat_lfc(params);
        act_design <- designMatrix(fit_options);
        act_contrast <- contrast(fit_options);
        act_adj_meth <- adjust_method(params);
        
        if (is(M, "IsoDataSet")) {
          isoC <- isoCounts(M);
          geneI <- isoGeneRel(M);
          expD <- expData(M);
          contrast <- levels(expD[,1]);
          expD[,1] <- as.factor(as.character(col_data(fit_options)[,]));
          
          colName <- colnames(expD)[[1]];
          isoDataSet <- IsoDataSet(isoC, expD, colName, geneI);
          isoDataSet <- buildLowExpIdx(isoDataSet);
          
          dsRes <- NBTest(isoDataSet, colName, test='F', contrast=contrast);
          res <- results(dsRes, filter=FALSE);
          resRank <- as.matrix(c(by(res, res$gene, function(x)
            sum(x$stat*sign(x$odd), na.rm=TRUE))));
          
          aux <- unique(res$gene);
          res <- res[, c('gene', 'genePval')];
          colnames(res) <- c('gene', 'p.value');
          res <- res[!is.na(res$p.value),];
          res <- unique(res);
          aux <- cbind(aux[!aux %in% res$gene], NA);
          colnames(aux) <- colnames(res);
          res <- rbind(res, aux);
          res <- res[match(rownames(M), res$gene),]; # reorder genes as in M
          # Adjusted pvalues
          res$p.adjust <- p.adjust(res$p.value, method=act_adj_meth);
          res$rank <- resRank[res$gene,];
        } else { # DGEList or MAList
          if (is(M, "DGEList")) { # apply voom
              M <- voom(M, design=act_design);
          }
          
          # Adjust the model
          fit <- lmFit(M, act_design);
          # treat correction
          fit2 <- treat(contrasts.fit(fit, act_contrast), lfc=act_treat_lfc);
          # Adjusted pvalues
          fit2$p.adjust <- apply(fit2$p.value, 2, p.adjust, method=act_adj_meth);
          res <- data.frame(gene=rownames(fit2), 
                            p.value=fit2$p.value, 
                            p.adjust=fit2$p.adjust,
                            rank=fit2$t);
        }
        rownames(res) <- res$gene;
        
        return(res);
    }
)

setGeneric(name="igsaGetDEGenes",
    def=function(seaParams, exprData, fitOptions) {
    standardGeneric("igsaGetDEGenes")
})

#'@importFrom futile.logger flog.info
#'@include ExprData-class.R
#'@include FitOptions-class.R
#'@include IGSAinput-class.R
setMethod(
    f="igsaGetDEGenes",
    signature=c("SEAparams", "ExprData", "FitOptions"),
    definition=function(seaParams, exprData, fitOptions) {
        # get the fit
        act_fit <- get_fit(exprData, fitOptions, seaParams);

        de_coff <- de_cutoff(seaParams);
        # get the DE genes depending on the cutoff value
        dif <- act_fit[,'p.adjust',drop=FALSE] <= de_coff;
        dif <- unique(rownames(dif)[dif]);
        
        flog.info(paste("DE genes", length(dif), "of a total of", 
                    nrow(exprData), "(", 
                    round(length(dif)/nrow(exprData)*100,2), "%)"));
        
        return(dif);
    }
)
