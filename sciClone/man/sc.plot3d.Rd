\name{sc.plot3d}
\alias{sc.plot3d}
\title{
  plot3d: Plot 
}
\description{
  Creates a 3d representation of clustering from a trio of input samples.
}
\usage{sc.plot3d(sco, samplesToPlot, size=700, outputFile=NULL)
}
\arguments{
  \item{sco}{
    an scObject object created by the sciClone() function
  }
  \item{samplesToPlot}{
    a character vector giving the names of the samples to plot. If only
    three samples were used in the clustering, just pass in sco@sampleNames
  }
  \item{size}{
    the size in pixels of the resulting square 3d plot
  }
  \item{outputFile}{
    if specified will create an animated gif that rotates the plot 360
    degrees. (Requires imagemagick to be installed on the host
    machine).
  }
}  
\details{  
  Uses rgl to create an interactive 3d plot of a trio of samples with
  clusters labelled.  Will also create an animated gif for display
  purposes (requires imagemagick to be installed on the host
  machine and 'convert' to be in the path).
}
\examples{
  # sc.plot3d(sc, samplesToPlot=c("sample1","sample2","sample3"), outputFile="results/clusters.3d.gif")
}
