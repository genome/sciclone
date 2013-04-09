\name{sc.plot2d}
\alias{sc.plot2d}
\title{
  plot2d: Plot 
}
\description{
  Creates a pdf with one plot for each pair of samples.
}
\usage{sc.plot2d(sco, outputFile, overlayClusters=TRUE, ellipse.metadata = list())
}
\arguments{
  \item{sco}{
    an scObject object created by the sciClone() function
  }
  \item{outputFile}{
    path/name of pdf file of plots to create
  }
  \item{overlayClusters}{
    boolean argument to decide whether to overlay cluster designations
    on the standard scatterplot output
  }
  \item{ellipse.metadata}{
    description of ellipses for plotting cluster std errors/deviations
  }
}  
\details{  
  Creates a scatter plot for each pairing of samples, with clusters
  indicated by color and shape.
}

\value{
  does not return anything, creates pdf file
}

\examples{
  # sc.plot2d(sc,"results/clusters.2d.pdf")
}
