setGeneric(name="IGSA", def=function(igsaInput, bp_param) {
    standardGeneric("IGSA")
})

#'@importClassesFrom BiocParallel BiocParallelParam
#'@importFrom futile.logger flog.info
#'@include DEnricher.R
#'@include IGSAinput.R
#'@include IGSAinput-class.R
#'@include IGSAinput-getterSetters.R
#'@include IGSAres.R
#'@include MGSZ.R
setMethod(
    f="IGSA",
    signature=c("IGSAinput", "BiocParallelParam"),
    definition=function(igsaInput, bp_param) {
        flog.info("*************************************");
        flog.info(paste(name(igsaInput), ": Starting IGSA analysis."));
        
        fit_options <- fitOptions(igsaInput);
        expr_data   <- exprData(igsaInput);
        
        actGeneSets <- geneSetsList(igsaInput);
        if (length(actGeneSets) == 0) {
            stop("No gene sets provided.");
        }
        
        # merging all gene sets in order to run just one mGSZ and one dEnricher,
        # its much faster, but after we will have to split these results
        merged_gene_sets <- merge_gene_sets(actGeneSets);
        flog.info(paste(length(geneSets(merged_gene_sets)), "Gene Sets."));
        
        flog.info(paste(name(igsaInput), ": dEnricher starting."));
        # run DEnricher
        deRes <- DEnricher(seaParams(igsaInput), expr_data, fit_options,
                        merged_gene_sets, useVoom(igsaInput), bp_param);
        flog.info(paste(name(igsaInput), ": dEnricher finnished."));
        
        flog.info(paste(name(igsaInput), ": mGSZ starting."));
        # run MGSZ
        mgszRes <- MGSZ(gseaParams(igsaInput), expr_data, fit_options,
                        merged_gene_sets, useVoom(igsaInput), bp_param);
        flog.info(paste(name(igsaInput), ": mGSZ finnished."));
        
        # splitting all results. it is a list of GenesetsRes objects
        splitted_res <- split_results(deRes, mgszRes, actGeneSets);
        
        # get the genes rank to give more information in the results object
        genes_rank <- get_fit(expr_data, fit_options, useVoom(igsaInput),
                                SEAparams());
        genes_rank <- data.frame(geneID=rownames(genes_rank),
                                    rank=genes_rank$t);
        colnames(genes_rank)[2] <- name(igsaInput);
        
        # create the IGSAres object
        igsaRes <- IGSAres(name=name(igsaInput), gene_sets_res=splitted_res,
                            genes_rank=genes_rank);
        
        flog.info(paste(name(igsaInput), ": IGSA analysis ended."));
        
        return(igsaRes);
    }
)

setGeneric(name="merge_gene_sets", def=function(geneSetsList) {
    standardGeneric("merge_gene_sets")
})

#'@include Genesets.R
#'@include Geneset.R
# input: list of Genesets
# output: a Genesets object
setMethod(
    f="merge_gene_sets",
    signature=c("list"),
    definition=function(geneSetsList) {
        stopifnot(all(unlist(lapply(geneSetsList, function(x)
                    is(x, "Genesets")))));
        
        # each individual gene set id will be its gene set name "_MIGSA_"
        # its id. e.g. BP gene set GO:10000 will be BP_MIGSA_GO:10000
        # Not the best idea haha
        merged <- Reduce(function(...) append(...),
            lapply(geneSetsList, function(gene_sets) {
                act_name <- name(gene_sets);
                act_gene_sets <-
                    lapply(geneSets(gene_sets), function(gene_set) {
                        id(gene_set) <- paste(act_name, "MIGSA", id(gene_set),
                                                sep="_");
                        return(gene_set);
                    }
                )
            })
        )
        return(Genesets(name="MIGSA_internal_use", gene_sets=merged));
    }
)

setGeneric(name="split_results",
    def=function(sea_res, gsea_res, geneSetsList) {
    standardGeneric("split_results")
})

#'@include GenesetRes.R
#'@include GenesetsRes.R
#'@include GSEAres.R
#'@include SEAres.R
# input geneSetsList is a list of Genesets
setMethod(
    f="split_results",
    signature=c("SEAres", "GSEAres", "list"),
    definition=function(sea_res, gsea_res, geneSetsList) {
        
        # So now lets split the merged results!!
        # and return to each its additional info they had (as is_GO)
        
        splitted_res <- lapply(geneSetsList, function(act_GS) {
            act_name <- name(act_GS);
            act_is_GO <- isGO(act_GS);
            # grep the ones which start with gene set name followed by "_MIGSA_"
            act_pattern <- paste("^", act_name, "_MIGSA_", sep="");
            act_sea_res <- lapply(gene_sets_res(sea_res), function(act_res) {
                if (grepl(act_pattern, id(act_res))) {
                    id(act_res) <- gsub(act_pattern, "", id(act_res));
                    return(act_res);
                } else {
                    return(NA);
                }
            })
            act_sea_res <- SEAres(
                            gene_sets_res=act_sea_res[!is.na(act_sea_res)]);
            
            act_gsea_res <- lapply(gene_sets_res(gsea_res), function(act_res) {
                if (grepl(act_pattern, id(act_res))) {
                    id(act_res) <- gsub(act_pattern, "", id(act_res));
                    return(act_res);
                } else {
                    return(NA);
                }
            })
            act_gsea_res <- GSEAres(
                            gene_sets_res=act_gsea_res[!is.na(act_gsea_res)]);
            
            return(GenesetsRes(gene_sets_name=act_name, is_GO=act_is_GO,
                                sea_res=act_sea_res, gsea_res=act_gsea_res));
        })
        return(splitted_res);
    }
)
