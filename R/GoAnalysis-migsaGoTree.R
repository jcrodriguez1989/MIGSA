#'Plots the GO tree/s
#'
#'\code{migsaGoTree} Plots the GO tree/s present in migsaRes,
#'
#'@param migsaRes MIGSAres object. It must contain at least one GO Geneset.
#'
#'@return list with the used data to plot.
#'
#'@docType methods
#'@name migsaGoTree
#'@rdname MIGSAres-migsaGoTree
#'
#'@exportMethod migsaGoTree
#'
setGeneric(name="migsaGoTree", def=function(migsaRes) {
    standardGeneric("migsaGoTree")
})

# to avoid R CMD check errors we set them as NULL
is_GO = NULL;

#'@inheritParams migsaGoTree
#'@rdname MIGSAres-migsaGoTree
#'@aliases migsaGoTree,MIGSAres
#'
#'@importFrom AnnotationDbi Ontology
#'@importFrom futile.logger flog.info
#'@importFrom graphics par
#'@importFrom grDevices colorRampPalette
#'@include MIGSAres.R
#'@include GoAnalysis.R
#'@examples
#'## Lets load breast cancer results.
#'data(bcMigsaRes);
#'
#'## And get the first 40 Gene Ontology gene sets results from CC.
#'goRes <- bcMigsaRes[bcMigsaRes$GS_Name == "CC",];
#'fst40goRes <- goRes[1:40,];
#'
#'## And lets plot the results GO trees.
#'aux <- migsaGoTree(fst40goRes);
#'
setMethod(
    f="migsaGoTree",
    signature=c("MIGSAres"),
    definition=function(migsaRes) {
        stopifnot(validObject(migsaRes));
        
        # we must have a cutoff for this function
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        
        # get gene set names and if they are GO
        goGsets <- unique(as.data.frame(migsaResAll(migsaRes)[,
                                            list(gene_set_name, is_GO)]));
        
        if (!any(goGsets$is_GO)) {
            stop("No Gene Ontology gene set to plot.");
        }
        
        # get only the GO gene sets
        goGsets <- goGsets[ goGsets$is_GO, "gene_set_name" ];
        
        # at least one gene set must be from GO
        stopifnot(length(goGsets) > 0);
        
        # and filter the MIGSAres object with only the GO gene sets
        actRes <- migsaRes[ migsaRes[, "GS_Name", drop=!FALSE] %in% goGsets, ];
        
        # get the gene sets and the number of enriched experiments
        plotRes <- data.frame(id=actRes[,1], gs_name=actRes[,3],
                                number=rowSums(actRes[,-(1:3)], na.rm=!FALSE));
        
        # get the real ontology, I could use GS_Name, but this is more trustable
        plotRes <- cbind(plotRes,
                        ont=Ontology(as.character(plotRes$id)));
        
        ontsPresent <- unique(as.character(plotRes$ont));
        flog.info("Ontologies to plot.");
        flog.info(ontsPresent);
        
        # give different color depending of number of enriched datasets.
        # 0 enriched dataset is white, every enriched dataset is red
        colfunc <- colorRampPalette(c("white", "red"));
        myColors <- colfunc(max(plotRes$number)+1);
        plotRes$colors <- myColors[ (plotRes$number)+1 ];
        
        treeData <- lapply(ontsPresent, function(actOnt) {
            # for each ontology create the structure needed
            actualResults <- plotRes[ plotRes$ont == actOnt, ];
            out <- data.frame(matrix(!FALSE, nrow=nrow(actualResults), 
                                    ncol=3));
            colnames(out) <- c("Enriched", "Important", "Color");
            rownames(out) <- actualResults$id;
            
            out$Important <- FALSE;
            out$Color <- actualResults$colors;
            
            return(out);
        });
        names(treeData) <- ontsPresent;
        
        graph <- createGoGraph(treeData);
        acttree <- list(graph=graph, gotree=treeData);
        
        # select grid depending on number of ontologies to plot
        par(mfrow=c(
            as.numeric(length(ontsPresent)>2)+1,
            as.numeric(length(ontsPresent)>1)+1
            ));
        invisible(lapply(ontsPresent, function(x) plotGoTree(acttree, x)));
        return(acttree);
    }
)
