##-------------------------------------------------
## Set up an object to hold the data,
##
initScClass <- function(){
  setClass("scObject", representation(clust="list", densities="list", dimensions="numeric",
                                      marginalClust="list", sampleNames="character", vafs.1d="list",
                                      vafs.merged="data.frame", ellipse.metadata="list", purities="numeric"))
}
