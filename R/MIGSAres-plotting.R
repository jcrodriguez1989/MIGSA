# to avoid R CMD check errors we set them as NULL
number = NULL;

#'MIGSAres plots
#'
#'\code{genesHeatmap} plots a heatmap with the number of experiments in which 
#'each gene contributed to enrich each gene set.
#'\code{genesBarplot} generates a barplot of the number of gene sets in which 
#'each gene contributed to enrich. Each gene set counts 1 regardless if it was 
#'enriched in many experiments. x-axis each gene, y-axis number of gene sets 
#'in which it contributed to enrich.
#'\code{migsaHeatmap} plots the enrichment heatmap of the MIGSAres object.
#'\code{geneSetBarplot} generates a barplot of the number of experiments in 
#'which each gene set was enriched. x-axis each gene set, y-axis times it was 
#'enriched (0 to #experiments).
#'
#'@param migsaRes MIGSAres object.
#'@param enrFilter numeric. Keep gene sets enriched in at least enrFilter 
#'experiments.
#'@param gsFilter numeric. Keep genes enriched in at least gsFilter gene sets.
#'@param expFilter numeric. Keep experiments which enriched at least expFilter 
#'gene sets.
#'@param col.dist character. Distance algorithm to be used in columns, passed 
#'to vegdist function.
#'@param row.dist character. Distance algorithm to be used in rows, passed to 
#'vegdist function.
#'@param ... In heatmap functions, ... are additional parameters passed to 
#'heatmap.2 function. In other functions it is not in use.
#'
#'@return In heatmap functions: A list returned by heatmap.2 function 
#'(plotted data). In other functions: A ggplot object used as graphic.
#'
#'@docType methods
#'@name MIGSAres-plots
#'@rdname MIGSAres-plots
#'@seealso \code{\link{setEnrCutoff}}
#'
#'@examples
#'data(migsaRes);
#'## For this example lets work with the first 50 results
#'migsaRes <- migsaRes[1:50,];
#'
#'#### genesHeatmap
#'## First lets set a cutoff of 0.1
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.1);
#'
#'## Lets check what genes contributed to enrich the highest number of gene 
#'## sets (in more than one gene set).
#'genesHeatmap(migsaResWCoff, gsFilter=1);
#'
#'## Moreover we can keep gene sets which where enriched in more than 
#'## enrFilter experiments. To do this, we can use the enrFilter parameter.
#'
#'#### genesBarplot
#'## Lets set a cutoff of 0.01
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.01);
#'
#'## Lets check what genes contributed to enrich the highest number of gene 
#'## sets (in more than one gene set).
#'genesBarplot(migsaResWCoff, gsFilter=1);
#'
#'## Moreover we can keep gene sets which where enriched in more than 
#'## enrFilter experiments. To do this, we can use the enrFilter parameter.
#'
#'#### migsaHeatmap
#'## Lets set a cutoff of 0.1
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.1);
#'
#'## Lets visually check enriched gene sets shared between experiments.
#'migsaHeatmap(migsaResWCoff);
#'
#'#### geneSetBarplot
#'## Lets set a cutoff of 0.1
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.1);
#'
#'## Lets check in how many experiments each gene set was enriched (in more 
#'## than one experiment).
#'geneSetBarplot(migsaResWCoff, enrFilter=1);
#'
setGeneric(name="MIGSAres-plots", def=function(migsaRes) {
    standardGeneric("MIGSAres-plots")
})

#'@name genesHeatmap
#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases genesHeatmap,MIGSAres-method
#'@exportMethod genesHeatmap
#'
setGeneric(name="genesHeatmap", def=function(migsaRes, ...) {
    standardGeneric("genesHeatmap")
})

#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases genesHeatmap,MIGSAres-method
#'
#'@importFrom gplots heatmap.2
#'@include MIGSAres-genesManipulation.R
#'@include MIGSAres-setEnrCutoff.R
#'
setMethod(f="genesHeatmap",
    signature=c("MIGSAres"),
    definition=function(migsaRes, enrFilter=0, gsFilter=0, ...) {
        stopifnot(validObject(migsaRes));
        
        # we must have a cutoff for this function
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        
        # keep gene sets enriched in more than enrFilter experiments
        actRes <- migsaRes[ rowSums(migsaRes[,-(1:3)], na.rm=!FALSE) >
                                    enrFilter, ];
        plotGenes <- genesInSets(actRes);
        
        # keep genes enriched in more than gsFilter gene sets
        plotGenes <- plotGenes[, colSums(plotGenes > 0) > gsFilter];
        
        p <- heatmap.2(plotGenes, labRow=rep("", nrow(plotGenes)), ...);
        p$data <- plotGenes;
        
        return(p);
    }
)

#'@name genesBarplot
#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases genesBarplot,MIGSAres-method
#'@exportMethod genesBarplot
#'
setGeneric(name="genesBarplot", def=function(migsaRes, ...) {
    standardGeneric("genesBarplot")
})

#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases genesBarplot,MIGSAres-method
#'
#'@importFrom ggplot2 aes geom_bar ggplot
#'@include MIGSAres-genesManipulation.R
#'@include MIGSAres-setEnrCutoff.R
#'
setMethod(
    f="genesBarplot",
    signature=c("MIGSAres"),
    definition=function(migsaRes, enrFilter=0, gsFilter=0) {
        stopifnot(validObject(migsaRes));
        
        # we must have a cutoff for this function
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        
        # keep gene sets enriched in more than enrFilter experiments
        actRes <- migsaRes[ rowSums(migsaRes[,-(1:3)], na.rm=!FALSE) >
                                    enrFilter, ];
        plotGenes <- genesInSets(actRes);
        
        # keep genes enriched in more than gsFilter gene sets
        plotGenes <- plotGenes[, colSums(plotGenes > 0) > gsFilter];
        
        # count 1 disregard if the gene enriched in more than 1 experiment
        plotGenes <- data.frame(id=colnames(plotGenes),
                                number=colSums(plotGenes > 0));
        
        ## todo: add symbol to plotGenes
        # plotGenes$symbol <- entrez2symbol(plotGenes$id);
        
        # reordering bars
        plotGenes$id <- factor(plotGenes$id,
                                levels=plotGenes$id[order(plotGenes$number,
                                                        decreasing=!FALSE)]);
        plotGenes <- plotGenes[order(plotGenes$number, decreasing=!FALSE),];
        
        p <- ggplot(plotGenes);
        p <- p + geom_bar(aes(id, y=number), stat="identity");
        print(p);
        
        return(p);
    }
)

#'@name migsaHeatmap
#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases migsaHeatmap,MIGSAres-method
#'@exportMethod migsaHeatmap
#'
setGeneric(name="migsaHeatmap", def=function(migsaRes, ...) {
    standardGeneric("migsaHeatmap")
})

#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases migsaHeatmap,MIGSAres-method
#'
#'@importFrom futile.logger flog.debug
#'@importFrom gplots heatmap.2
#'@importFrom grDevices rainbow
#'@importFrom stats as.dendrogram hclust
#'@importFrom vegan vegdist
#'@include MIGSAres.R
#'@include MIGSAres-setEnrCutoff.R
#'
setMethod(
    f="migsaHeatmap",
    signature=c("MIGSAres"),
    definition=function(migsaRes, enrFilter=0, expFilter=0,
    col.dist="jaccard", row.dist=col.dist, ... ) {
        otherParams <- list(...);
        # we must have a cutoff for this function
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        
        allRes <- get_summary(migsaRes);
        
        # experiment names (present in migsaRes)
        exp_cols <- setdiff(colnames(allRes), c("id", "Name", "GS_Name"));
        
        # terms filtering by enrFilter and expFilter
        keepRows <- rowSums(allRes[,exp_cols], na.rm=!FALSE)/
                        ncol(allRes[,exp_cols]) >= enrFilter;
        keepCols <- colSums(allRes[,exp_cols], na.rm=!FALSE)/
                        nrow(allRes[,exp_cols]) >= expFilter;
        allRes <- allRes[ keepRows, keepCols ];
        
        # however rows with no enrichment by any method are deleted
        allRes <- allRes[rowSums(allRes[,exp_cols], na.rm=!FALSE) > 0,];
        flog.debug(paste("In migsaHeatmap, after filtering, dim=",
                        dim(allRes)));
        
        # get color for each different gene set
        rowColors <- rep("white", nrow(allRes));
        if ("GS_Name" %in% colnames(allRes)) {
            rowColors <- as.character(allRes$GS_Name);
            pal <- rainbow(length(unique(rowColors)));
            names(pal) <- unique(rowColors);
            rowColors <- pal[rowColors]; rm(pal);
        }
        if ("RowSideColors" %in% names(otherParams)) {
            rowColors <- otherParams[["RowSideColors"]];
            otherParams <- otherParams[names(otherParams) != "RowSideColors"];
#           rowColors <- RowSideColors; ## ojo! ver si funca
        }
        
        numRes <- apply(allRes[,exp_cols], 2, as.numeric);
        
        colColors <- rep("white", length(exp_cols));
#         if (!is.null(conditions)) {
#             pal <- rainbow(length(conditions));
#             for (i in 1:length(conditions)) {
#                 colColors[ grep(conditions[[i]], exp_cols) ] <- pal[[i]];
#                 colnames(numRes) <- gsub(conditions[[i]], "",
#                         colnames(numRes));
#             }
#         }
        if ("ColSideColors" %in% names(otherParams)) {
            colColors <- otherParams[["ColSideColors"]];
            otherParams <- otherParams[names(otherParams) != "ColSideColors"];
#           colColors <- ColSideColors; ## ojo! ver si funca
        }
        
        # jaccard clustering per column (experiment)
        dd.s <- suppressWarnings(vegdist(t(numRes), col.dist, na.rm=!FALSE));
        dd.s[is.na(dd.s)] <- 0;
        h.s <- hclust(dd.s, method="average");
        
        # jaccard clustering per row (gene set)
        ddr.s <- suppressWarnings(vegdist(numRes, row.dist, na.rm=!FALSE));
        ddr.s[is.na(ddr.s)] <- 0;
        hr.s <- hclust(ddr.s, method="average");
        
        # we use values -1 NA (not analyzed), 0 not enriched, 1 enriched
        numRes[is.na(numRes)] <- -1;
        allParams <- c(otherParams, list(x=as.matrix(numRes), 
            Rowv=as.dendrogram(hr.s), labRow=rep("", nrow(numRes)), 
            Colv=as.dendrogram(h.s), colsep=1:(ncol(numRes)-1), 
            sepwidth=c(0.025, 0.025), RowSideColors=rowColors, 
            ColSideColors=colColors));
        p <- do.call(heatmap.2, allParams);
        p$data <- numRes;
        # how to add legends:
        # http://stackoverflow.com/questions/17041246/
        #    how-to-add-an-inset-subplot-to-topright-of-an-r-plot
        
        return(p);
    }
)

#'@name geneSetBarplot
#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases geneSetBarplot,MIGSAres-method
#'@exportMethod geneSetBarplot
#'
setGeneric(name="geneSetBarplot", def=function(migsaRes, ...) {
    standardGeneric("geneSetBarplot")
})

#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases geneSetBarplot,MIGSAres-method
#'
#'@importFrom ggplot2 aes geom_bar ggplot
#'@include MIGSAres-setEnrCutoff.R
#'
setMethod(
    f="geneSetBarplot",
    signature=c("MIGSAres"),
    definition=function(migsaRes, enrFilter=0) {
        stopifnot(validObject(migsaRes));
        
        # we must have a cutoff for this function
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        
        # keep gene sets enriched in more than enrFilter experiments
        actRes <- migsaRes[ rowSums(migsaRes[,-(1:3)], na.rm=!FALSE) >
                                enrFilter, ];
        plotRes <- data.frame(actRes[,1:3],
                                number=rowSums(actRes[,-(1:3)], na.rm=!FALSE));
        
        # reordering bars
        plotRes$id <- factor(plotRes$id,
                            levels=plotRes$id[order(plotRes$number,
                                                    decreasing=!FALSE)]);
        plotRes <- plotRes[order(plotRes$number, decreasing=!FALSE),];
        
        p <- ggplot(plotRes);
        p <- p + geom_bar(aes(id, y=number), stat="identity");
        print(p);
        
        return(p);
    }
)
