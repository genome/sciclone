##-------------------------------------------------
## Set up an object to hold the data,
##
initScClass <- function(){
  setClass("scObject", representation(clust="list", densities="list", dimensions="numeric",
                                      marginalClust="list", sampleNames="character", vafs.1d="list",
                                      vafs.merged="data.frame", ellipse.metadata="list"))
  ## setMethod("initialize", "scObject",
  ##         function(.Object){
  ##           .Object@clust=as.list(c())
  ##           .Object@densities=list(c())
  ##           .Object@dimensions=numeric()
  ##           .Object@readInfo=data.frame()
  ##           .Object@chrs=as.list(c())
  ##           return(.Object)              
  ##         })
}
