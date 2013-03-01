\name{sciClone}
\alias{sciClone}
\title{
  sciClone: Illustrates sub-clones within a sequenced sample
}
\description{
  sciClone integrates the read depth and copy number information at single nucleotide variant locations in order to create a visualization of the sub-clonal variant populations areas of deleted, neutral, and amplified copy number. Variants are optionally clustered to formalize predictions of the sub-clonal makeup of the sample.
}
\usage{
    sciClone(vafs, outputImage, copyNumberCalls = NULL, sampleName, 
            minimumDepth = 100, clusteredDataOutputFile = NULL, componentDistribution = "Binomial", 
            maximumClustersToTest = 8, minimumLabelledPeakHeight = 0.001, 
            onlyLabelHighestPeak = FALSE, overlayClusters = FALSE, plotOnlyCN2 = FALSE, 
            positionsToHighlight = NULL, purity = 0, highlightSexChrs = TRUE, 
            testing = FALSE, cnCallsAreLog2 = FALSE, useSexChrs = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vafs}{
    variant allele frequency data for single nucleotide variants in
    5-column format: 1) chromosome 2) position 3) reference supporting
    read counts 4) variant supporting read counts 5) variant allele
    frequency
  }
  \item{outputImage}{
    filename of the PDF output image
  }
  \item{copyNumberCalls}{
    windowed copy number designations in 4-column format: 1) chromosome
    2) window start position 3) window stop position 4) copy number
    designation within window; unrepresented regions are assumed to have
    copy number equal to 2
  }
  \item{regionsToExclude}{
    windowed exclusion regions in 3-column format: 1) chromosome 2) window
    start position 3) window stop position; single nucleotide variants
    falling into these windows will not be included in the analysis; use
    for copy-number-neutral LOH regions, for example
  }
  \item{sampleName}{
    text describing the data or sample (used for plotting purposes)
  }
  \item{minimumDepth}{
    read depth cutoff used for excluding lower-read-depth variants
  }
  \item{clusteredDataOutputFile}{
    output file containing cluster designations for clustered data
    points
  }
  \item{componentDistribution}{
    distribution used in modeling variant subclones (options are "Beta",
    "Normal", or "Binomial")
  }
  \item{maximumClustersToTest}{
    max number of clusters to consider when chosing the best component
    fit to the data
  }
  \item{minimumLabelledPeakHeight}{
    minimum density at which to label a peak in the output kernel
    density plot
  }
  \item{onlyLabelHighestPeak}{
    boolean argument to decide whether to label only the highes peak in
    the output kernel density plot
  }
  \item{overlayClusters}{
    boolean argument to decide whether to overlay cluster designations
    on the standard scatterplot output
  }
  \item{plotOnlyCN2}{
    boolean argument to decide weather to plot a variant scatterplot for
    all copy number designations (FALS), or only the neutral points
    (TRUE)
  }
  \item{positionsToHighlight}{
    a list of positions in 2-column format 1) chromosome 1) position to
    highlight on the output scatterplot by using a different plotting
    emblem than that used for other data points (a black x is currently
    used)
  }
  \item{purity}{
    user-provided purity estimation of a tumor sample (integer from 1 to
    100)
  }
  \item{highlightSexChrs}{
    boolean argument to decide whether or not to highlight all
    sex-chromosome variants using a different plotting emblem than that
    used for other data points (a black x is currently used)
  }
  \item{testing}{
    boolean argument specifying the preference to additionally output
    all component fits in separate images
  }
  \item{cnCallsAreLog2}{
    boolean argument specifying whether or not the copy number
    predictions are in log2 format (as opposed to being absolute copy
    number designations)
  }
  \item{useSexChrs}{
    boolean argument to specify preference of whether (TRUE) or not
    (FALSE) to use variants on sex chromosomes in the plotting and
    clustering steps of the tool.
  }
  \item{doClustering}{
    boolean argument - if (TRUE), the tool will attempt to use clustering
    to identify subclones. If (FALSE) this stage is skipped, and only
    the plots are produced.
  }
  \item{highlightsHaveNames}{
    boolean argument - if (TRUE), the tool will look for a third column
    in the file specified by positionsToHighlight and will number these 
    points then create a legend with the names. Also sets plotOnlyCN2 to
    true, to make room for the legend.
  }

}
\details{
  Testing details.
}

\value{
  Testing value.
}

\examples{

}
