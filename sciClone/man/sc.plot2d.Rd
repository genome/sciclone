\name{sc.plot2d}
\alias{sc.plot2d}
\title{
  plot2d: Plot 
}
\description{
  Creates a pdf with one plot for each pair of samples.
}
\usage{sc.plot2d(sco, outputFile, positionsToHighlight=NULL,
                 highlightsHaveNames=FALSE, overlayClusters=TRUE,
                 overlayErrorBars=FALSE, ellipse.metadata = list(),
                 singlePage=FALSE, scale=1)
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
  \item{overlayErrorBars}{
    boolean argument to decide whether to overlay (binomial) error bars
    on VAFs
  }
  \item{positionsToHighlight}{
    a list of positions in 2-column format 1) chromosome 2) position.
    These points will be highlighted on the output scatterplot by using
    a different character than that used for other data points
  }
  \item{highlightsHaveNames}{
    boolean argument - if (TRUE), the tool will look for a third column
    in the file specified by positionsToHighlight and will label these
    points on the plot.
  }
  \item{ellipse.metadata}{
    description of ellipses for plotting cluster std errors/deviations
  }
  \item{singlePage}{
    If true, condenses all output plots onto a single page.
  }
  \item{scale}{
    scale the spacing of the plots and fonts up or down to produce output
    suitable for presentations (try 1.0) or print (try 2.0).
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
