#'Accessors for MIGSAinput class.
#'
#'Getters and setters functions to access MIGSAinput object slots.
#'
#'@param object MIGSAinput object.
#'@param value value to replace in the slot.
#'
#'@return modified object or desired slot.
#'
#'@docType methods
#'@name MIGSAinput-getterSetters
#'@rdname MIGSAinput-getterSetters
#'@include MIGSAinput-class.R
#'
#'@examples
#'## Lets create an empty MIGSAinput object and access its slots.
#'migsaInput <- MIGSAinput();
#'
#'## Lets get out migsaInput values.
#'bpParam(migsaInput);
#'experiments(migsaInput);
#'geneSetsList(migsaInput);
#'
setGeneric(name="MIGSAinput-getterSetters", def=function(object) {
    standardGeneric("MIGSAinput-getterSetters")
})

#'@name experiments
#'@rdname MIGSAinput-getterSetters
#'@inheritParams MIGSAinput-getterSetters
#'@aliases experiments,MIGSAinput-method
#'@exportMethod experiments
#'
setGeneric(name="experiments", def=function(object) {
    standardGeneric("experiments")
})

#'@rdname MIGSAinput-getterSetters
#'@inheritParams MIGSAinput-getterSetters
#'@aliases experiments,MIGSAinput-method
#'
setMethod(f="experiments", signature="MIGSAinput",
    definition=function(object) {
        return(object@experiments)
    }
)

#'@name experiments<-
#'@rdname MIGSAinput-getterSetters
#'@inheritParams MIGSAinput-getterSetters
#'@aliases experiments<-,MIGSAinput-method
#'@exportMethod experiments<-
#'
setGeneric(name="experiments<-", def=function(object, value) {
    standardGeneric("experiments<-")
})

#'@rdname MIGSAinput-getterSetters
#'@inheritParams MIGSAinput-getterSetters
#'@aliases experiments<-,MIGSAinput-method
#'
setReplaceMethod(f="experiments", signature="MIGSAinput",
    definition=function(object, value) {
        object@experiments <- value;
        validObject(object);
        return(object);
    }
)

#'@name bpParam
#'@rdname MIGSAinput-getterSetters
#'@inheritParams MIGSAinput-getterSetters
#'@aliases bpParam,MIGSAinput-method
#'@exportMethod bpParam
#'
setGeneric(name="bpParam", def=function(object) {
    standardGeneric("bpParam")
})

#'@rdname MIGSAinput-getterSetters
#'@inheritParams MIGSAinput-getterSetters
#'@aliases bpParam,MIGSAinput-method
#'
setMethod(f="bpParam", signature="MIGSAinput",
    definition=function(object) {
        return(object@bp_param)
    }
)

#'@name bpParam<-
#'@rdname MIGSAinput-getterSetters
#'@inheritParams MIGSAinput-getterSetters
#'@aliases bpParam<-,MIGSAinput-method
#'@exportMethod bpParam<-
#'
setGeneric(name="bpParam<-", def=function(object, value) {
    standardGeneric("bpParam<-")
})

#'@rdname MIGSAinput-getterSetters
#'@inheritParams MIGSAinput-getterSetters
#'@aliases bpParam<-,MIGSAinput-method
#'
setReplaceMethod(f="bpParam", signature="MIGSAinput",
    definition=function(object, value) {
        object@bp_param <- value;
        validObject(object);
        return(object);
    }
)

#'@rdname MIGSAinput-getterSetters
#'@inheritParams MIGSAinput-getterSetters
#'@aliases geneSetsList,MIGSAinput-method
#'@include IGSAinput.R
#'
setMethod(f="geneSetsList", signature="MIGSAinput",
    definition=function(object) {
        return(object@gene_sets_list)
    }
)

#'@rdname MIGSAinput-getterSetters
#'@inheritParams MIGSAinput-getterSetters
#'@aliases geneSetsList<-,MIGSAinput-method
#'
setReplaceMethod(f="geneSetsList", signature="MIGSAinput",
    definition=function(object, value) {
        object@gene_sets_list <- value;
        validObject(object);
        return(object);
    }
)
