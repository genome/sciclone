\name{writeClusterTable}
\alias{writeClusterTable}
\title{
  writeClusterTable
}
\description{
  creates a tab-delimted file of variants and clusters
}
\usage{writeClusterTable(sco, outputFile)
}
\arguments{
  \item{sco}{
    an scObject object created by the sciClone() function
  }
  \item{outputFile}{
    path/name of tab-delimted file to create
  }
}  
\details{  
  creates a tab-delimted file that combines variants, allele
  frequencies, cluster assignments and probabilities.
}
\examples{
  # writeClusterTable(sc, "/path/to/outputFile.tsv")
}
