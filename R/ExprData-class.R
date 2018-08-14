#'ExprData S4 class implementation in R
#' 
#'This S4 class contains the expression data, it can be a MAList, DGEList,
#'or IsoDataSet.
#'Important: Rownames are going to be used as gene identifiers, make sure that 
#'this IDs are the same used in your gene sets.
#'If x is an IsoDataSet, then gene_id column from isoGeneRel(x) will be used
#'as rownames.
#'If the ExprData is a matrix, then MIGSA will use it only for mGSZ. And it 
#'will consider first column as the gene's ranking for the real data, and each
#'remaining column for the permutations.
#'
#'@docType methods
#'@name ExprData-class
#'@rdname ExprData-class
#'
#'@importClassesFrom edgeR DGEList
#'@importClassesFrom limma MAList
#'@importClassesFrom NBSplice IsoDataSet
#'@exportClass ExprData
#'@examples
#'## Lets create a ExprData object from counts data (DGEList).
#'rnaData <- matrix(rnbinom(10000,mu=5,size=2),ncol=4);
#'rownames(rnaData) <- 1:nrow(rnaData); # It must have rownames (gene names).
#'
#'## Now we can use rnaExprData as an ExprData object.
#'\dontrun{
#'rnaExprData <- DGEList(counts=rnaData); 
#'}
#'
#'## Lets create a ExprData object from micro array data (MAList).
#'maData <- matrix(rnorm(10000),ncol=4);
#'rownames(maData) <- 1:nrow(maData); # It must have rownames (gene names).
#'
#'## Now we can use maExprData as an ExprData object.
#'maExprData <- new("MAList",list(M=maData));
#'
setClassUnion("ExprData", c("MAList", "DGEList", "IsoDataSet", "matrix"));
