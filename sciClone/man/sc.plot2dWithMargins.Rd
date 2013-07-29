\name{sc.plot2dWithMargins}
\alias{sc.plot2dWithMargins}
\title{
  plot2d: Plot
}
\description{
  Creates a pdf where each page has a scatterplot for each pair of
  samples along with a pair of adjacent histograms for each sample.
}
\usage{sc.plot2dWithMargins(sco, outputFile, positionsToHighlight=NULL,
highlightsHaveNames=FALSE, overlayErrorBars=FALSE)
}
\arguments{
  \item{sco}{
    an scObject object created by the sciClone() function
  }
  \item{outputFile}{
    path/name of pdf file of plots to create
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
  \item{overlayErrorBars}{
    boolean argument to decide whether to overlay (binomial) error bars
    on VAFs
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
