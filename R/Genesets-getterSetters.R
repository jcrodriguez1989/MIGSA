#'Accessors for Genesets class.
#'
#'Getters and setters functions to access Genesets object slots.
#'
#'@param object Genesets object.
#'@param value value to replace in the slot.
#'
#'@return modified object or desired slot.
#'
#'@docType methods
#'@name Genesets-getterSetters
#'@rdname Genesets-getterSetters
#'@include Genesets-class.R
#'
#'@examples
#'## Lets create manually a gene sets containing three gene set.
#'## First lets create manually three gene set.
#'myGs1 <- Geneset(id="fakeId1", name="fakeName1", genes=as.character(1:10));
#'myGs2 <- Geneset(id="fakeId2", name="fakeName2", genes=as.character(7:15));
#'myGs3 <- Geneset(id="fakeId3", name="fakeName2", genes=as.character(20:28));
#'
#'## Now we can create our Genesets object.
#'myGSs <- Genesets(name="myGenesets", gene_sets=list(myGs1, myGs2, myGs3),
#'is_GO=FALSE);
#'
#'## Lets get out Genesets values, and modify its name.
#'name(myGSs);
#'name(myGSs) <- "newName";
#'geneSets(myGSs);
#'isGO(myGSs);
#'
setGeneric(name="Genesets-getterSetters", def=function(object) {
    standardGeneric("Genesets-getterSetters")
})

#'@name name
#'@inheritParams Genesets-getterSetters
#'@rdname Genesets-getterSetters
#'@aliases name,Genesets-method
#'@exportMethod name
#'
setGeneric(name="name", def=function(object) {
    standardGeneric("name")
})

#'@rdname Genesets-getterSetters
#'@inheritParams Genesets-getterSetters
#'@aliases name,Genesets-method
#'
setMethod(f="name", signature="Genesets",
    definition=function(object) {
        return(object@name)
    }
)

#'@name name<-
#'@rdname Genesets-getterSetters
#'@inheritParams Genesets-getterSetters
#'@aliases name<-,Genesets-method
#'@exportMethod name<-
#'
setGeneric(name="name<-", def=function(object, value) {
    standardGeneric("name<-")
})

#'@rdname Genesets-getterSetters
#'@inheritParams Genesets-getterSetters
#'@aliases name<-,Genesets-method
#'
setReplaceMethod(f="name", signature="Genesets",
    definition=function(object, value) {
        object@name <- value;
        validObject(object);
        return(object);
    }
)

#'@name geneSets
#'@inheritParams Genesets-getterSetters
#'@rdname Genesets-getterSetters
#'@aliases geneSets,Genesets-method
#'@exportMethod geneSets
#'
setGeneric(name="geneSets", def=function(object) {
    standardGeneric("geneSets")
})

#'@rdname Genesets-getterSetters
#'@inheritParams Genesets-getterSetters
#'@aliases geneSets,Genesets-method
#'
setMethod(f="geneSets", signature="Genesets",
    definition=function(object) {
        return(object@gene_sets)
    }
)

#'@name geneSets<-
#'@inheritParams Genesets-getterSetters
#'@rdname Genesets-getterSetters
#'@aliases geneSets<-,Genesets-method
#'@exportMethod geneSets<-
#'
setGeneric(name="geneSets<-", def=function(object, value) {
    standardGeneric("geneSets<-")
})

#'@inheritParams Genesets-getterSetters
#'@rdname Genesets-getterSetters
#'@aliases geneSets<-,Genesets-method
#'
setReplaceMethod(f="geneSets", signature="Genesets",
    definition=function(object, value) {
        object@gene_sets <- value
        validObject(object)
        return(object)
    }
)

#'@name isGO
#'@inheritParams Genesets-getterSetters
#'@rdname Genesets-getterSetters
#'@aliases isGO,Genesets-method
#'@exportMethod isGO
#'
setGeneric(name="isGO", def=function(object) {
    standardGeneric("isGO")
})

#'@rdname Genesets-getterSetters
#'@inheritParams Genesets-getterSetters
#'@aliases isGO,Genesets-method
#'
setMethod(f="isGO", signature="Genesets",
    definition=function(object) {
        return(object@is_GO)
    }
)

#'@name isGO<-
#'@inheritParams Genesets-getterSetters
#'@rdname Genesets-getterSetters
#'@aliases isGO<-,Genesets-method
#'@exportMethod isGO<-
#'
setGeneric(name="isGO<-", def=function(object, value) {
    standardGeneric("isGO<-")
})

#'@inheritParams Genesets-getterSetters
#'@rdname Genesets-getterSetters
#'@aliases isGO<-,Genesets-method
#'
setReplaceMethod(f="isGO", signature="Genesets",
    definition=function(object, value) {
        object@is_GO <- value
        validObject(object)
        return(object)
    }
)
