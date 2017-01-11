setGeneric(name="DEnricher", def=function(params, M, fit_options, gene_sets,
                                            use_voom, bp_param) {
    standardGeneric("DEnricher")
})

#'@importClassesFrom BiocParallel BiocParallelParam
#'@importFrom futile.logger flog.info
#'@include FitOptions.R
#'@include Genesets.R
#'@include Geneset.R
#'@include IGSAinput.R
#'@include SEAparams.R
#'@include SEAres.R
# params <- seaParams(igsaInput); M <- expr_data; 
# use_voom <- useVoom(igsaInput); gene_sets <- merged_gene_sets
setMethod(
    f="DEnricher",
    signature=c("SEAparams", "ExprData", "FitOptions", "Genesets", "logical",
                                                "BiocParallelParam"),
    definition=function(params, M, fit_options, gene_sets, use_voom, 
                        bp_param) {
        if (length(de_genes(params)) == 1 && is.na(de_genes(params)[[1]])) {
            dif <- igsaGetDEGenes(params, M, fit_options, use_voom);
        } else {
            dif <- de_genes(params);
            flog.info(paste("DE genes", length(dif), "of a total of", nrow(M),
                            "(", round(length(dif)/nrow(M)*100,2), "%)"));
        }
        
        br <- br(params);
        allGenes <- unique(unlist(asList(gene_sets)));
        
        if (length(br) > 1) {
            br <- intersect(br, allGenes);
            flog.info(paste("Using user provided BR:", length(br), "genes."));
            if (length(br) < 2) {
                stop("No genes in br after intersecting with experiment genes");
            }
        } else if (br == "bri") {
            br <- allGenes;
            flog.info(paste("Using BRI:", length(br), "genes."));
        } else if (br == "briii") {
            br <- unique(rownames(M));
            flog.info(paste("Using BRIII:", length(br), "genes."));
        } else {
            stop("Incorrect br option.");
        }
        
        dif <- intersect(dif, br);
        validGSets <- lapply(geneSets(gene_sets), function(actGs) {
            validGenes <- intersect(genes(actGs), br);
            if (length(validGenes) == 0) {
                return(NULL);
            }
            genes(actGs) <- validGenes;
            return(actGs);
        });
        validGSets <- validGSets[ !unlist(lapply(validGSets, is.null)) ];
        geneSets(gene_sets) <- validGSets;
        
        test <- test(params);
        seaRes <- runDEnricher(dif, gene_sets, br, test=test, bp_param);
        
        seaRes <- SEAres(gene_sets_res=seaRes);
        
        return(seaRes);
    }
)

setGeneric(name="runDEnricher",
    def=function(dif, genesets, br, test, bp_param) {
    standardGeneric("runDEnricher")
})

#'@importClassesFrom BiocParallel BiocParallelParam
#'@importFrom BiocParallel bplapply
#'@importFrom futile.logger flog.info
#'@importFrom stats fisher.test pbinom phyper
#'@include GenesetRes.R
#'@include Geneset.R
#'@include Genesets.R
# genesets <- gene_sets; test="FisherTest"
setMethod(
    f="runDEnricher",
    signature=c("character", "Genesets", "character", "character",
                                                    "BiocParallelParam"),
    definition=function(dif, genesets, br, test, bp_param) {
        GeneID <- dif;
        genes.group <- GeneID[!is.na(GeneID)];
        
        gs <- genesets;
        nSet <- length(geneSets(gs))
        
        doFisherTest <- function(genes.group, genes.term, genes.universe) {
            genes.hit <- intersect(genes.group, genes.term)
            X <- length(genes.hit)
            K <- length(genes.group)
            M <- length(genes.term)
            N <- length(genes.universe)
            cTab <- matrix(c(X, K - X, M - X, N - M - K + X), nrow = 2, 
                dimnames = list(c("anno", "notAnno"),
                                c("group", "notGroup")))
            p.value <- ifelse(all(cTab == 0), 1, stats::fisher.test(cTab, 
                                            alternative = "greater")$p.value)
            return(p.value)
        }
        doHypergeoTest <- function(genes.group, genes.term, genes.universe) {
            genes.hit <- intersect(genes.group, genes.term)
            X <- length(genes.hit)
            K <- length(genes.group)
            M <- length(genes.term)
            N <- length(genes.universe)
            x <- X
            m <- M
            n <- N - M
            k <- K
            p.value <- ifelse(m == 0 || k == 0, 1, stats::phyper(x, 
                                    m, n, k, lower.tail=FALSE, log.p=FALSE))
            return(p.value)
        }
        doBinomialTest <- function(genes.group, genes.term, genes.universe) {
            genes.hit <- intersect(genes.group, genes.term)
            X <- length(genes.hit)
            K <- length(genes.group)
            M <- length(genes.term)
            N <- length(genes.universe)
            p.value <- ifelse(K == 0 || M == 0 || N == 0, 1, stats::pbinom(X, 
                                    K, M/N, lower.tail=FALSE, log.p=FALSE))
            return(p.value)
        }
        
        genes.universe <- br;
        genes.group <- intersect(genes.universe, genes.group)
        if (length(genes.group) == 0) {
            warning(paste("There are no differentialy expressed genes",
                "being used."));
        }
        
        flog.info(paste("Running SEA at cores:", bp_param$workers));
        seaRes <- bplapply(geneSets(gs), function(actGset) {
        #     seaRes <- lapply(gs@gene_sets, function(actGset) {
            genes.term <- unique(unlist(genes(actGset)))
            p.value <- switch(test,
                FisherTest   = doFisherTest(genes.group, genes.term,
                                                            genes.universe),
                HypergeoTest = doHypergeoTest(genes.group, genes.term,
                                                            genes.universe),
                BinomialTest = doBinomialTest(genes.group, genes.term,
                                                            genes.universe))
            
            impGenes <- intersect(genes.group, genes.term);
            actRes <- GenesetRes(id=id(actGset), name=getName(actGset),
                                pvalue=p.value, genes=genes(actGset),
                                enriching_genes=impGenes);
            return(actRes);
        }, BPPARAM=bp_param)
        #     })

        return(seaRes);
    }
)
