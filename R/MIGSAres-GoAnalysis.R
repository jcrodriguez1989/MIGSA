# to avoid R CMD check errors we set them as NULL
is_GO = NULL;

#'Explore Gene Ontology gene sets in MIGSAres
#'
#'\code{migsaGoTree} plots the GO tree/s present in migsaRes.
#'\code{getHeights} returns the heights of given a list of ids (GO IDs).
#'
#'@param migsaRes MIGSAres object. It must contain at least one GO gene set.
#'@param ids character vector indicating the queried GO ids.
#'@param minHeight logical indicating if the minimum or maximum height must be 
#'calculated. If it is FALSE then the longest path to the root is calculated, 
#'otherwise, the shortest path.
#'@param ... not in use.
#'
#'@return If migsaGoTree: A list with the used data to plot. If getHeights: A 
#'list with each term height.
#'
#'@docType methods
#'@name MIGSAres-GOanalysis
#'@rdname MIGSAres-GoAnalysis
#'
#'@examples
#'## Lets load breast cancer results.
#'data(bcMigsaRes);
#'
#'###### migsaGoTree
#'## Get the first 40 Gene Ontology gene sets results from CC.
#'goRes <- bcMigsaRes[bcMigsaRes$GS_Name == "CC",];
#'fst40goRes <- goRes[1:40,];
#'
#'## And lets plot the results GO trees.
#'\dontrun{
#'aux <- migsaGoTree(fst40goRes);
#'}
#'
#'###### getHeights
#'## Get the first 40 Gene Ontology gene sets IDs.
#'goIds <- bcMigsaRes[bcMigsaRes$GS_Name %in% c("BP", "CC", "MF"), "id"];
#'fst40goIds <- goIds[1:40,];
#'
#'\dontrun{
#'## And lets get the heights in the GO tree structure.
#'getHeights(fst40goIds);
#'}
#'
setGeneric(name="MIGSAres-GOanalysis", def=function(migsaRes) {
    standardGeneric("MIGSAres-GOanalysis")
})

#'@name migsaGoTree
#'@inheritParams MIGSAres-GOanalysis
#'@rdname MIGSAres-GoAnalysis
#'@aliases migsaGoTree,MIGSAres-method
#'@exportMethod migsaGoTree
#'
setGeneric(name="migsaGoTree", def=function(migsaRes) {
    standardGeneric("migsaGoTree")
})


#'@inheritParams MIGSAres-GOanalysis
#'@rdname MIGSAres-GoAnalysis
#'@aliases migsaGoTree,MIGSAres-method
#'
#'@importFrom AnnotationDbi Ontology
#'@importFrom futile.logger flog.info
#'@importFrom graphics par
#'@importFrom grDevices colorRampPalette
#'@include MIGSAres.R
#'@include GoAnalysis.R
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

#'@name getHeights
#'@inheritParams MIGSAres-GOanalysis
#'@rdname MIGSAres-GoAnalysis
#'@aliases getHeights,character-method
#'@exportMethod getHeights
#'
setGeneric(name="getHeights", def=function(ids, ...) {
    standardGeneric("getHeights")
})

#'@inheritParams MIGSAres-GOanalysis
#'@rdname MIGSAres-GoAnalysis
#'@aliases getHeights,character-method
#'
#'@importFrom GO.db GOBPPARENTS GOCCPARENTS GOMFPARENTS
#'@include GoAnalysis.R
#'
setMethod(
    f="getHeights",
    signature=c("character"),
    definition=function(ids, minHeight=TRUE) {
        allParents <- as.list(GOBPPARENTS);
        allParents <- c(allParents, as.list(GOMFPARENTS));
        allParents <- c(allParents, as.list(GOCCPARENTS));
        
        ids <- gsub(" ", "", ids);
        
        # for each GO id find its height
        result <- lapply(ids, function(x) {
            getHeight(x, minHeight, allParents)
        });
        
        return(unlist(result));
    }
)
