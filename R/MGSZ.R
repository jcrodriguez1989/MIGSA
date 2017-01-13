setGeneric(name="MGSZ", def=function(params, M, fit_options, genesets,
                                    use_voom, bp_param) {
    standardGeneric("MGSZ")
})

# params <- gseaParams(igsaInput); M <- exprData(igsaInput);
# genesets <- merged_gene_sets; use_voom <- useVoom(igsaInput)
#'@importClassesFrom BiocParallel BiocParallelParam
#'@include FitOptions.R
#'@include Geneset.R
#'@include GenesetRes.R
#'@include Genesets.R
#'@include GSEAparams.R
#'@include GSEAres.R
#'@include IGSAinput.R
setMethod(
    f="MGSZ",
    signature=c("GSEAparams", "ExprData", "FitOptions", "Genesets", "logical",
                "BiocParallelParam"),
    definition=function(params, M, fit_options, genesets, use_voom,
                        bp_param) {
        # mGSZ uses the gene sets as a list
        mgsz.gene.sets <- asList(genesets);
        
        # check the ranking function depending if we use voom or not
        if(use_voom) {
            rankFunc <- voomLimaRank;
        } else {
            rankFunc <- mGszEbayes;
        }
        
        # run faster version of mGSZ
        mgszRes <- MIGSA_mGSZ(M, fit_options, mgsz.gene.sets, rankFunc, params,
                                bp_param);
        
        # convert mGSZ results (data.frame) to GenesetRes
        mgszRes <- lapply(geneSets(genesets), function(actGset) {
            mgszGSres <- mgszRes[ mgszRes$gene.sets == id(actGset), ];
            if (nrow(mgszGSres) == 0) {
                pval <- as.numeric(NA);
                actScore <- as.numeric(NA);
                actImpGenes <- "";
            } else {
                pval <- mgszGSres$pvalue;
                actScore <- mgszGSres$mGszScore;
                actImpGenes <- as.character(mgszGSres$impGenes);
                # todo: put better genes separator
                actImpGenes <- strsplit(actImpGenes, "_sep_")[[1]];
            }
            
            actRes <- GenesetRes(id=id(actGset), name=getName(actGset),
                                score=actScore, pvalue=pval,
                                genes=genes(actGset),
                                enriching_genes=actImpGenes);
        })
        
        # this is a GSEAres
        gseaRes <- GSEAres(gene_sets_res=mgszRes);
        
        return(gseaRes);
    }
)

setGeneric(name="MIGSA_mGSZ",
    def=function(M, fit_options, gene.sets, rankFunction, params,
                        bp_param) {
    standardGeneric("MIGSA_mGSZ")
})

# gene.sets <- mgsz.gene.sets; rankFunction <- rankFunc
# M <- igsaInput@expr_data; fit_options <- igsaInput@fit_options;
#'@importClassesFrom BiocParallel BiocParallelParam
#'@importClassesFrom edgeR DGEList
#'@importFrom BiocParallel bplapply
#'@importFrom data.table data.table
#'@importFrom edgeR [.DGEList
#'@importFrom futile.logger flog.debug flog.info
#'@importFrom mGSZ count_hyge_var_mean count.prob.sum mGSZ.p.values 
#'rm.small.genesets toMatrix
#'@include FitOptions.R
#'@include GSEAparams.R
#'@include IGSAinput.R
setMethod(
    f="MIGSA_mGSZ",
    signature=c("ExprData", "FitOptions", "list", "function", "GSEAparams",
                                            "BiocParallelParam"),
    definition=function(M, fit_options, gene.sets, rankFunction, params,
                                                    bp_param) {
        if (is(M, "DGEList")) {
            stopifnot(all(rownames(M$samples) == colnames(M$counts)));
        }
        
        ## all this code is from mGSZ, with some improvements
        
        expr.data <- M;
        min.cl.sz <- minSz(params); # min gene set size
        pre.var <- pv(params);
        wgt1 <- w1(params);
        wgt2 <- w2(params);
        var.constant <- vc(params);
        start.val=5;
        perm.number <- perm_number(params);
        
        flog.debug(paste("mGSZ: perm.number=", perm.number));
        flog.debug(paste("mGSZ: gene.sets number=", length(gene.sets)));
        
        expr.dat.sz <- dim(expr.data)
        # convert the gene sets list to a data.frame (genes x genesets)
        gene.sets <- toMatrix(expr.data, gene.sets)
        num.genes <- nrow(gene.sets)
        if(!(num.genes == expr.dat.sz[1])){
            stop("Number of genes in gene set data and expression data
                    do not match")
        }
        rm(expr.dat.sz);
        
        # Remove genes which are not in any term
        data <- rowSums(gene.sets) > 0;
        expr.data <- expr.data[data, ];
        gene.sets <- gene.sets[data, ];
        rm(data);
        
        # Remove too small gene sets
        aux <- ncol(gene.sets);
        gene.sets <- rm.small.genesets(gene.sets, min.cl.sz)
        flog.debug(paste("mGSZ: Not analyzed gene sets:",
                    aux-ncol(gene.sets)));
        rm(aux);
        
        # now convert it to data.table to make it much faster
        gene.sets <- data.table(gene.sets)
        
        if (ncol(expr.data) < 1 || nrow(expr.data) < 1 || 
            ncol(gene.sets) < 1 || nrow(gene.sets) < 1) {
            stop("No gene sets to test after mGSZ filtering");
        }
        
        geneset.classes <- colnames(gene.sets)
        num.classes <- ncol(gene.sets)
        num.genes <- nrow(expr.data)
        geneset.class.sizes <- colSums(gene.sets)
        
        # Calculation of gene scores
        tmp = MIGSA_diffScore(expr.data, fit_options, perm.number, 
                        rankFunction, bp_param);
        
        tmp.expr.data = tmp$rankings
        perm.number = tmp$perm.number
        diff.expr.dat.sz <- dim(tmp.expr.data)
        
        # calculate hyge and prob (mGSZ stuff)
        set_sz <- geneset.class.sizes;
        hyge_stat <- count_hyge_var_mean(nrow(gene.sets), unique(set_sz));
        prob_sum <- lapply(1:length(unique(set_sz)), function(j) {
            count.prob.sum(nrow(gene.sets), hyge_stat$mean[, j],
                hyge_stat$var[, j])
        });
        
        # mGSZ score for real data
        pos.scores <- MIGSA_mGSZ.test.score2(tmp.expr.data[,1], gene.sets,
                                            rownames(expr.data), wgt1, wgt2,
                                            pre.var, var.constant, start.val,
                                            set_sz, hyge_stat, prob_sum);
        
        pos.mGSZ.scores <- as.numeric(pos.scores$mgszScore);
        
        flog.info(paste("Running score at cores:", bp_param$workers));
        
        # mGSZ score for permuted data
        col.perm.mGSZ <- do.call(rbind, 
            bplapply(1:(diff.expr.dat.sz[2] -1), function(k) {
        #   lapply(1:(diff.expr.dat.sz[2] -1), function(k) {
                tmp <- MIGSA_mGSZ.test.score(tmp.expr.data[, k+1], gene.sets,
                                        wgt1, wgt2, pre.var, var.constant,
                                        start.val, set_sz, hyge_stat,
                                        prob_sum);
                flog.debug(paste("Finnished perm n:", k));
                return(tmp);
        #   })
            }, BPPARAM=bp_param)
        );
        
        # p-value calculation for the gene set scores
        mGSZ.p.vals.col.perm <- mGSZ.p.values(pos.mGSZ.scores,col.perm.mGSZ)
        
        # preparing output table
        mGSZ.table.col <- data.frame(
            gene.sets=colnames(gene.sets)[mGSZ.p.vals.col.perm$class.ind],
            set.size=geneset.class.sizes[mGSZ.p.vals.col.perm$class.ind],
            gene.set.scores=pos.mGSZ.scores[mGSZ.p.vals.col.perm$class.ind],
            pvalue=mGSZ.p.vals.col.perm$EV.class,
            mGszScore=pos.scores[mGSZ.p.vals.col.perm$class.ind, "rawMgsz"],
            impGenes=pos.scores[mGSZ.p.vals.col.perm$class.ind, "impGenes"],
            row.names=NULL
        )
        
        # sorting the results table
        out <- mGSZ.table.col[order(mGSZ.table.col$pvalue,decreasing=FALSE),]
        
        return(out)
    }
)

setGeneric(name="MIGSA_diffScore",
    def=function(data, fit_options, perm.number, rankFunction,
                bp_param) {
    standardGeneric("MIGSA_diffScore")
})

#'@importClassesFrom BiocParallel BiocParallelParam
#'@importFrom BiocParallel bplapply
#'@importFrom futile.logger flog.info
#'@include IGSAinput.R
#'@include FitOptions.R
setMethod(
    f="MIGSA_diffScore",
    signature=c("ExprData", "FitOptions", "numeric", "function",
                "BiocParallelParam"),
    definition=function(data, fit_options, perm.number, rankFunction,
                        bp_param) {
        # data: expression data matrix
        # labels: Vector of response values (example: 1,2)
        # perms.number: Number of sample permutations
        
        dime2 <- dim(data)
        pit <- ncol(data)
        
        # generate permutations
        all_perms <- replicate(perm.number, sample(1:pit, replace=FALSE));
        
        # it deletes duplicate permutations
        unique.perm <- unique(t(all_perms))
        if (dim(unique.perm)[1] < perm.number) {
            flog.info(paste("Number of unique permutations:",
                dim(unique.perm)[1]));
            all_perms <- t(unique.perm)
        }
        dime  <- dim(all_perms)
        dime2 <- dim(data)
        rankings = array(0, c(dime2[1], dime[2] + 1))
        
        # gene scores for real data
        rankings[,1] <- rankFunction(data, fit_options);
        
        flog.info(paste("Getting ranking at cores:", bp_param$workers));
        
        # gene scores for permuted data
        calcRes <- bplapply(1:ncol(all_perms), function(i) {
        #     calcRes <- lapply(1:ncol(all_perms), function(i) {
            # modify the design matrix using the order given by the permutation
            permDesign <- designMatrix(fit_options)[all_perms[,i],];
            new_fit_options <- fit_options;
            designMatrix(new_fit_options) <- permDesign;
            
            actRank <- rankFunction(data, new_fit_options);
            
            # todo: maybe the best alternative would be delete the gene from
            # everywhere
            if (any(is.na(actRank))) {
                warning(paste(sum(is.na(actRank)), 
                    "genes generated NAs when estimating fit for perm", i));
                actRank[is.na(actRank)] <- mean(actRank, na.rm=!FALSE);
                # genes which are NA are given the mean value of all genes, 
                # so they dont contribute to enrichment score.
            }
            
            return(actRank);
        }, BPPARAM=bp_param);
        #     })
        calcRes <- matrix(unlist(calcRes), ncol=ncol(all_perms));
        
        rankings[,-1] <- calcRes;
        
        out <- list(rankings=rankings, perm.number=dim(unique.perm)[1])
        return(out);
    }
)

#'@importFrom limma eBayes contrasts.fit lmFit
#'@include FitOptions.R
mGszEbayes <- function(exprMatrix, fit_options) {
#     flog.info("Using ebayes");
    design <- designMatrix(fit_options);
    contrast <- contrast(fit_options);
    fit1 <- lmFit(exprMatrix, design);
    fit2 <- eBayes(contrasts.fit(fit1, contrast));
    
    res <- fit2$t;
    return(res);
}

## voom + limma
#'@importFrom limma eBayes contrasts.fit lmFit voom
#'@include FitOptions.R
voomLimaRank <- function(exprMatrix, fit_options) {
#     flog.info("Using voom+limma");
    design <- designMatrix(fit_options);
    contrast <- contrast(fit_options);
    newExpr <- voom(exprMatrix, design);
    
    # Adjust the model
    fit1 <- lmFit(newExpr, design);
    fit2 <- eBayes(contrasts.fit(fit1, contrast));
    
    res <- fit2$t;
    return(res);
}

setGeneric(name="MIGSA_mGSZ.test.score2",
    def=function(expr.data, gene.sets, genes, wgt1, wgt2, pre.var, 
                var.constant, start.val, set_sz, hyge_stat, prob_sum) {
    standardGeneric("MIGSA_mGSZ.test.score2")
})

# expr.data <- tmp.expr.data[,1]; genes <- rownames(expr.data); 
#'@importFrom data.table data.table
setMethod(
    f="MIGSA_mGSZ.test.score2",
    signature=c("numeric", "data.table", "character", "numeric", "numeric",
                "numeric", "numeric", "numeric", "numeric", "list", "list"),
    definition=function(expr.data, gene.sets, genes, wgt1, wgt2, pre.var,
                    var.constant, start.val, set_sz, hyge_stat, prob_sum) {
        num.genes <- length(expr.data)
        # order the data by genes decreasing value
        ord_out <- order(expr.data, decreasing=TRUE)
        expr.data <- expr.data[ord_out]
        genes <- genes[ord_out];
        gene.sets <- gene.sets[ord_out,] 
        
        expr.data.ud <- expr.data[num.genes:1]
        num.classes=ncol(gene.sets)
        unique_class_sz_ln <- length(unique(set_sz))
        
        # calculating some mGSZ stuff
        pre_z_var.1 <- MIGSA_sumVarMean_calc(expr.data, gene.sets, pre.var,
                                            set_sz, hyge_stat, prob_sum)
        pre_z_var.2 <- MIGSA_sumVarMean_calc(expr.data.ud, gene.sets, pre.var,
                                            set_sz, hyge_stat, prob_sum)
        
        # calculating some mGSZ stuff
        Z_var1 = MIGSA_calc_z_var(num.genes, unique_class_sz_ln,
                                pre_z_var.1$Z_var, wgt2, var.constant)
        Z_var2 = MIGSA_calc_z_var(num.genes, unique_class_sz_ln,
                                pre_z_var.2$Z_var, wgt2, var.constant)
        
        # for each gene set calculate enrichment score
        out <- do.call(rbind, lapply(1:num.classes, function(k) {
            # genes in set and out set
            po1 <- which(gene.sets[[k]] == 1)
            po0 <- which(gene.sets[[k]] == 0)
            
            tmp1 = expr.data
            tmp1[po0] = 0
            tmp0 = expr.data
            tmp0[po1] = 0
            
            result1 = cumsum(tmp1 - tmp0) -
            pre_z_var.1$Z_mean[, pre_z_var.1$class_size_index[k]]
            result2 = cumsum(tmp1[num.genes:1] - tmp0[num.genes:1]) - 
            pre_z_var.2$Z_mean[, pre_z_var.2$class_size_index[k]]
            
            result1[1:start.val] <- 0
            result2[1:start.val] <- 0
            
            # A genes are in the same order as expr.data
            # b genes are backwards
            A = result1/Z_var1[, pre_z_var.1$class_size_index[k]]
            B = result2/Z_var2[, pre_z_var.2$class_size_index[k]]
            
            eScores <- c(A,B);
            maxPoint <- which.max(abs(c(A,B)));
            rawMgsz <- eScores[[maxPoint]];
            
            # important genes are the ones until the top of the enrichment 
            # score. If ES is positive then the fst ones, if negative, the last 
            # ones
            if (maxPoint <= length(A)) {
                impGenes <- intersect(genes[po1], genes[1:maxPoint]);
            } else {
                impGenes <- intersect(genes[po1],
                                rev(genes)[1:(maxPoint%%length(A))]);
            }
            
            # todo: define a better gene separator
            impGenes <- paste(impGenes, collapse="_sep_");
            res <- data.frame(mgszScore=abs(rawMgsz), rawMgsz=rawMgsz,
                                impGenes=impGenes);
            
            return(res);
        }));
        
        return(out)
    }
)

setGeneric(name="MIGSA_mGSZ.test.score",
    def=function(expr.data, gene.sets, wgt1, wgt2, pre.var, var.constant,
                start.val, set_sz, hyge_stat, prob_sum) {
    standardGeneric("MIGSA_mGSZ.test.score")
})

#'@importFrom data.table data.table
setMethod(
    f="MIGSA_mGSZ.test.score",
    signature=c("numeric", "data.table", "numeric", "numeric", "numeric",
                "numeric", "numeric", "numeric", "list", "list"),
    definition=function(expr.data, gene.sets, wgt1, wgt2, pre.var, 
            var.constant, start.val, set_sz, hyge_stat, prob_sum) {
        num.genes <- length(expr.data)
        # order the data by genes decreasing value
        ord_out <- order(expr.data, decreasing=TRUE)
        expr.data <- expr.data[ord_out]
        gene.sets <- gene.sets[ord_out, ]
        
        expr.data.ud <- expr.data[num.genes:1]
        num.classes = ncol(gene.sets)
        unique_class_sz_ln <- length(unique(set_sz))
        
        # calculating some mGSZ stuff
        pre_z_var.1 <- MIGSA_sumVarMean_calc(expr.data, gene.sets, pre.var,
                                        set_sz, hyge_stat, prob_sum)
        pre_z_var.2 <- MIGSA_sumVarMean_calc(expr.data.ud, gene.sets, pre.var,
                                        set_sz, hyge_stat, prob_sum)
        
        # calculating some mGSZ stuff
        Z_var1 = MIGSA_calc_z_var(num.genes, unique_class_sz_ln,
                                pre_z_var.1$Z_var, wgt2, var.constant)
        Z_var2 = MIGSA_calc_z_var(num.genes, unique_class_sz_ln,
                                pre_z_var.2$Z_var, wgt2, var.constant)
        
        # for each gene set calculate enrichment score
        out <- unlist(lapply(1:num.classes, function(k) {
            po1 <- which(gene.sets[[k]] == 1)
            po0 <- which(gene.sets[[k]] == 0)
            
            #         if (length(po1) > 0 & length(po0) > 0) {
            tmp1 = expr.data
            tmp1[po0] = 0
            tmp0 = expr.data
            tmp0[po1] = 0
            
            result1 = cumsum(tmp1 - tmp0) -
                pre_z_var.1$Z_mean[, pre_z_var.1$class_size_index[k]]
            result2 = cumsum(tmp1[num.genes:1] - tmp0[num.genes:1]) - 
                pre_z_var.2$Z_mean[, pre_z_var.2$class_size_index[k]]
            
            result1[1:start.val] <- 0
            result2[1:start.val] <- 0
            
            # A genes are in the same order as expr.data
            # b genes are backwards
            A = result1/Z_var1[, pre_z_var.1$class_size_index[k]]
            B = result2/Z_var2[, pre_z_var.2$class_size_index[k]]
            
            return(max(abs(c(A, B))));
        }));
        
        return(out)
    }
)

setGeneric(name="MIGSA_sumVarMean_calc",
    def=function(expr_data, gene.sets, pre.var, set_sz, hyge_stat, prob_sum) {
    standardGeneric("MIGSA_sumVarMean_calc")
})

setMethod(
    f="MIGSA_sumVarMean_calc",
    signature=c("numeric", "data.frame", "numeric", "numeric", "list", "list"),
    definition=function(expr_data, gene.sets, pre.var, set_sz, hyge_stat,
                        prob_sum) {
        # this is a mGSZ function, I just put lapply instead of for loop
        dim_sets <- dim(gene.sets)
        unique_class_sz <- unique(set_sz)
        num_genes <- length(expr_data)
        divider <- 1:num_genes
        mean_table <- cumsum(expr_data)/divider
        mean_table_sq <- mean_table^2
        var_table <- cumsum(expr_data^2)/divider - (mean_table)^2 + pre.var
        
        class_sz_index <- do.call(c, lapply(1:dim_sets[2], function(i) {
            which(set_sz[i] == unique_class_sz)
        }));
        
        var_table <- var_table + pre.var
        
        max_value <- 1:dim_sets[1]
        calcRes <- lapply(1:length(unique_class_sz), function(j) {
            act_prob_sum <- prob_sum[[j]];
            z_var_j <- 4 * (var_table * act_prob_sum + mean_table_sq * 
                            hyge_stat$var[, j])
            z_mean_j <- mean_table * (2 * hyge_stat$mean[, j] - max_value)
            return(list(var=z_var_j, mean=z_mean_j));
        })
        
        z_var  <- do.call(cbind, lapply(calcRes, function(x) x$var));
        z_mean <- do.call(cbind, lapply(calcRes, function(x) x$mean));
        
        out = list(Z_var=z_var, Z_mean=z_mean, 
            class_size_index=class_sz_index, var_table=var_table, 
            mean_table_sq=mean_table_sq, set_sz=set_sz)
        return(out)
    }
)

setGeneric(name="MIGSA_calc_z_var",
    def=function(num.genes, unique_class_sz_ln, pre_z_var, wgt2,
            var.constant) {
    standardGeneric("MIGSA_calc_z_var")
})

#'@importFrom matrixStats colMedians
setMethod(
    f="MIGSA_calc_z_var",
    signature=c("integer", "integer", "matrix", "numeric", "numeric"),
    definition=function(num.genes, unique_class_sz_ln, pre_z_var, wgt2,
                        var.constant) {
        # this is a mGSZ function, I just put colMedias instead of for loop
        ones_matrix = matrix(1, num.genes, unique_class_sz_ln)
        ones_tmatrix = t(ones_matrix)
        median_matrix = colMedians(pre_z_var)
        pre_median_part = ones_tmatrix * median_matrix * wgt2
        median_part = t(pre_median_part)
        z_var = (pre_z_var + median_part + var.constant)^0.5
        return(z_var);
    }
)
