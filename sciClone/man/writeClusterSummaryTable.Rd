\name{writeClusterSummaryTable}
\alias{writeClusterSummaryTable}
\title{
  writeClusterSummaryTable
}
\description{
  Creates a table that lists each cluster along with it's mean VAF, and
  the standard error of the mean (high and low)
}
\usage{writeClusterSummaryTable(sco, outputFile)
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
}
\examples{
  # sc.plot3d(sc, samplesToPlot=c("sample1","sample2","sample3"), outputFile="results/clusters.3d.gif")
}
