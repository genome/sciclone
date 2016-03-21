\name{sc.plot2d}
\alias{sc.plot2d}
\title{
  plot2d: Plot 
}
\description{
  Creates a pdf with one plot for each pair of samples.
}
\usage{sc.plot2d(sco, outputFile=NULL, positionsToHighlight=NULL,
                 highlightsHaveNames=FALSE, overlayClusters=TRUE,
                 overlayErrorBars=FALSE, ellipse.metadata = list(),
                 singlePage=FALSE, scale=1, xlim=100, ylim=100,
                 plot.title=NULL, samplesToPlot=NULL,
                 clusterLegend=TRUE, 
                 flipSamples=FALSE, xlab=NULL,
                 ylab=NULL, colors=NULL,
                 plotWidth=7.2, plotHeight=6)

}
\arguments{
  \item{sco}{
    an scObject object created by the sciClone() function
  }
  \item{outputFile}{
    path/name of pdf file that will contain the plots. Optional, but 
    generally recommended, especially if you have more than two samples
  }
  \item{samplesToPlot}{
    a list of sample names that should be included in the 2d
    plots. Default is to plot allcombinations of all samples.
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
  \item{xlim}{
    max value of the xaxis
  }
  \item{ylim}{
    max value of the yaxis
  }
  \item{plot.title}{
    Set a new title that will be applied to all plots. The default is to 
    add a title with format "sampleNameA vs sampleNameB" to each plot
  }
  \item{xlab}{
    label for the x axis
  } 
  \item{ylab}{
    label for the y axis
  }
  \item{colors}{
    provide your own list of colors to use for the clusters
  }
  \item{flipSamples}{
    flip the samples being plotted on each axis (x moves to y and vice-versa)
  }
 \item{clusterLegend}{
   add a legend showing cluster numbers (TRUE/FALSE)
 }
 \item{plotWidth}{
   width of the pdf (in inches) that is produced (default 7.2)
 }
 \item{plotHeight}{
   height of the pdf (in inches)that is produced (default 6)
 }
  
}  
\details{  
  Creates a scatter plot for each pairing of samples, with clusters
  indicated by color and shape.
}

\value{
  does not return anything, creates pdf file or output plot
}

\examples{
  # sc.plot2d(sc,"results/clusters.2d.pdf")
}
