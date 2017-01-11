#'Plots enrichment heatmap
#'
#'\code{migsaHeatmap} plots the enrichment heatmap of the MIGSAres object.
#'
#'@param migsaRes MIGSAres object.
#'@param flevel numeric. Keeps gene sets enriched in at least flevel 
#'experiments.
#'@param penrich numeric. Keeps experiments which enriched at least penrich 
#'gene sets.
#'@param col.dist character. Distance algorithm to be used in columns, passed 
#'to vegdist function.
#'@param row.dist character. Distance algorithm to be used in rows, passed to 
#'vegdist function.
#'@param ... additional parameters passed to heatmap.2 function.
#'
#'@return list returned by heatmap.2 function (plotted data).
#'
#'@docType methods
#'@name migsaHeatmap
#'@rdname MIGSAres-migsaHeatmap
#'
#'@include MIGSAres-class.R
#'@exportMethod migsaHeatmap
#'
setGeneric(name="migsaHeatmap", def=function(migsaRes, ...) {
    standardGeneric("migsaHeatmap")
})

#'@inheritParams migsaHeatmap
#'@rdname MIGSAres-migsaHeatmap
#'@aliases migsaHeatmap,MIGSAres
#'
#'@importFrom futile.logger flog.debug
#'@importFrom gplots heatmap.2
#'@importFrom grDevices rainbow
#'@importFrom stats as.dendrogram hclust
#'@importFrom vegan vegdist
#'@include MIGSAres.R
#'@include MIGSAres-setEnrCutoff.R
#'@seealso \code{\link{setEnrCutoff}}
#'@examples
#'data(migsaRes);
#'
#'## First lets set a cutoff of 0.1
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.1);
#'
#'## Lets visually check enriched gene sets shared between experiments.
#'migsaHeatmap(migsaResWCoff);
#'
setMethod(
    f="migsaHeatmap",
    signature=c("MIGSAres"),
    definition=function(migsaRes, flevel=0, penrich=0,
    col.dist="jaccard", row.dist=col.dist, ... ) {
        otherParams <- list(...);
        
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        
        allRes <- get_summary(migsaRes);
        exp_cols <- setdiff(colnames(allRes), c("id", "Name", "GS_Name"));
        
        # terms filtering by flevel and penrich
        keepRows <- rowSums(allRes[,exp_cols], na.rm=!FALSE)/
                        ncol(allRes[,exp_cols]) >= flevel;
        keepCols <- colSums(allRes[,exp_cols], na.rm=!FALSE)/
                        nrow(allRes[,exp_cols]) >= penrich;
        allRes <- allRes[ keepRows, keepCols ];
        
        # however rows with no enrichment by any method are deleted
        allRes <- allRes[rowSums(allRes[,exp_cols], na.rm=!FALSE) > 0,];
        flog.debug(paste("In migsaHeatmap, after filtering, dim=",
                        dim(allRes)));
        
        # get color for each different database
        rowColors <- rep("white", nrow(allRes));
        if ("GS_Name" %in% colnames(allRes)) {
            rowColors <- as.character(allRes$GS_Name);
            pal <- rainbow(length(unique(rowColors)));
            names(pal) <- unique(rowColors);
            rowColors <- pal[rowColors]; rm(pal);
        }
        if ("RowSideColors" %in% names(otherParams)) {
            rowColors <- otherParams[["RowSideColors"]];
#           rowColors <- RowSideColors; ## ojo! ver si funca
        }
        
        # jaccard clustering per column (database)
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
#           colColors <- ColSideColors; ## ojo! ver si funca
        }
        
        dd.s <- suppressWarnings(vegdist(t(numRes), col.dist, na.rm=!FALSE));
        dd.s[is.na(dd.s)] <- 0;
        h.s <- hclust(dd.s, method="average");
        
        # jaccard clustering per row (term)
        ddr.s <- suppressWarnings(vegdist(numRes, row.dist, na.rm=!FALSE));
        ddr.s[is.na(ddr.s)] <- 0;
        hr.s <- hclust(ddr.s, method="average");
        
        numRes[is.na(numRes)] <- -1;
        p <- heatmap.2(as.matrix(numRes), Rowv=as.dendrogram(hr.s),
                    labRow=rep("", nrow(numRes)), Colv=as.dendrogram(h.s),
                    colsep=1:(ncol(numRes)-1), sepwidth=c(0.025, 0.025),
                    RowSideColors=rowColors, ColSideColors=colColors, ...);
        p$data <- numRes;
        # ver como agregar legends
        # http://stackoverflow.com/questions/17041246/
        #    how-to-add-an-inset-subplot-to-topright-of-an-r-plot
        
        return(p);
    }
)
