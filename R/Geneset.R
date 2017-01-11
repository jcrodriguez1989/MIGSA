#'@include Geneset-class.R
setGeneric(name="id", def=function(object) {
    standardGeneric("id")
})

setMethod(f="id", signature="Geneset",
    definition=function(object) {
        return(object@id)
    }
)

setGeneric(name="id<-", def=function(object, value) {
    standardGeneric("id<-")
})

setReplaceMethod(f="id", signature="Geneset",
    definition=function(object, value) {
        object@id <- value
        validObject(object)
        return(object)
    }
)

# backslah because it is internal
setGeneric(name="getName", def=function(object) {
    standardGeneric("getName")
})

setMethod(f="getName", signature="Geneset",
    definition=function(object) {
        return(object@name)
    }
)

setGeneric(name="genes", def=function(object) {
    standardGeneric("genes")
})

setMethod(f="genes", signature="Geneset",
    definition=function(object) {
        return(object@genes)
    }
)

setGeneric(name="genes<-", def=function(object, value) {
    standardGeneric("genes<-")
})

setReplaceMethod(f="genes", signature="Geneset",
    definition=function(object, value) {
        object@genes <- value
        validObject(object)
        return(object)
    }
)
