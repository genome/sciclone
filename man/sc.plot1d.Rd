\name{sc.plot1d}
\alias{sc.plot1d}
\title{
  plot1d: Plot kernel density and copy number information for each sample
}
\description{
  Creates a pdf with one plot for each sample.
}
\usage{sc.plot1d(sco, outputFile=NULL, cnToPlot=c(1,2,3,4),
                showCopyNumberScatterPlots=TRUE, highlightSexChrs=TRUE,
                positionsToHighlight=NULL, highlightsHaveNames=FALSE,
                overlayClusters=TRUE, overlayIndividualModels=TRUE,
                showHistogram=FALSE, showTitle=TRUE,
                biggerText=FALSE, highlightsOnHistogram=FALSE)

}
\arguments{
  \item{sco}{
    an scObject object created by the sciClone() function
  }
  \item{outputFile}{
    path/name of plot-containing pdf file to create
  }
  \item{cnToPlot}{
    create a shortened plot with only the density plots and indicated cn
    regions (a subset of c(1,2,3,4))
  }
  \item{showCopyNumberScatterPlots}{
    show the scatter plots of depth and VAF
  }
  \item{highlightSexChrs}{
    boolean argument to decide whether or not to highlight all
    sex-chromosome variants using a different plotting emblem than that
    used for other data points (a triangle is currently used)
  }
  \item{positionsToHighlight}{
    a list of positions in 2-column format 1) chromosome 2) position.
    These points will be highlighted on the output scatterplot by using
    a different character than that used for other data points
  }
  \item{highlightsHaveNames}{
    boolean argument - if (TRUE), the tool will look for a third column
    in the file specified by positionsToHighlight and will number these
    points, then create a legend with the names. Also sets cnToPlot=2
    to make room for the legend.
  }
  % \item{minimumLabelledPeakHeight}{
  %   minimum density at which to label a peak in the output kernel
  %   density plot
  % }
  % \item{onlyLabelHighestPeak}{
  %   boolean argument to decide whether to label only the highest peak in
  %   the output kernel density plot
  % }
  \item{overlayClusters}{
    boolean argument to decide whether to overlay cluster designations
    on the standard scatterplot output
  }
  \item{overlayIndividualModels}{
    show the model fits for each cluster on the density plot
  }
  \item{showHistogram}{
    display a histogram on the density plot
  }
  \item{showTitle}{
    Use the sample names as the titles of the plots
  }
  \item{biggerText}{
    bump up the font size to make small plots legible
  }
  \item{highlightsOnHistogram}{
    show higlighted points on histogram
  }
}
\value{
  does not return anything, creates pdf file
}
\examples{
  # sc.plot1d(sc,"results/clusters.1d.pdf")
}
