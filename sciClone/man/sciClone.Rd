\name{sciClone}
\alias{sciClone}
\title{
  sciClone: Identifies sub-clones within a sequenced sample
}
\description{
  sciClone integrates the read depth and copy number information at
  single nucleotide variant locations and clusters the variants in
  copy neutral regions, to formalize description of the sub-clonal
  architecture of the sample.
}
\usage{
    sciClone(vafs, copyNumberCalls=NULL, regionsToExclude=NULL,
             sampleNames, minimumDepth=100, clusterMethod="bmm", clusterParams=NULL,
             purities=NULL, cnCallsAreLog2=FALSE, useSexChrs=TRUE,
             doClustering=TRUE, verbose=TRUE, copyNumberMargins=0.25,
             maximumClusters=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{vafs}{
    a list of dataframes containing variant allele frequency data for
    single nucleotide variants in 5-column format:
    1. chromosome 2. position 3. reference-supporting
    read counts 4. variant-supporting read counts 5. variant allele
    frequency (between 0-100)
  }
  
  \item{copyNumberCalls}{
    list of dataframes containing windowed copy number segments in
    4-column format: 1. chromosome  2. window start position 3. window
    stop position 4. copy number value for that segment. Unrepresented
    regions are assumed to have a copy number of 2.
  }
  \item{regionsToExclude}{
    windowed exclusion regions in 3-column format: 1. chromosome 2. window
    start position 3. window stop position; Single nucleotide variants
    falling into these windows will not be included in the analysis. Use
    this input for LOH regions, for example
  }
  \item{sampleNames}{
    vector of names describing the samples
  }
  \item{minimumDepth}{
    cutoff used for excluding low-depth variants
  }
  \item{maximumClusters}{
    max number of clusters to consider when chosing the best component
    fit to the data.
  }
  \item{purities}{
    Purity estimation for a tumor sample (integer from 1 to 100). If not
    provided, the tool will attempt to infer this value from the data.
  }
  \item{cnCallsAreLog2}{
    boolean argument specifying whether or not the copy number
    predictions are in log2 format (as opposed to being absolute copy
    number designations)
  }
  \item{useSexChrs}{
    boolean argument to specify preference of whether (TRUE) or not
    (FALSE) to use variants on sex chromosomes in the clustering
    steps of the tool.
  }
  \item{doClustering}{
    boolean argument - if (TRUE), the tool will attempt to use clustering
    to identify subclones. If (FALSE) this stage is skipped, and an
    object suitable for feeding into the plotting functions is produced.    
  }
  \item{clusterMethod}{
    The framework is in place to drop in different clustering methods,
    but only 'bmm' is currently available.
  }
  \item{clusterParams}{
    The framework is in place to drop in different clustering methods
    and provide them with additional parameters, but only 'bmm' is
    currently available, so this should stay NULL.
  }
  \item{verbose}{
    if TRUE, prints lots of output to the screen that might be useful
    for debugging.
  }
  \item{copyNumberMargins}{
    In order to identify cleanly copy-number neutral regions, sciClone
    only considers sites with a copy number of 2.0 +/- this value. For
    example, if set to 0.25, regions at 2.20 will be considered cn-neutral,
    and regions at, 2.35 will not.
  }
}
\details{
  Testing details.
}

\value{
  returns a sciClone object containing merged vafs, clusters, and other
  information needed for visualization
}

\examples{
  #sc = sciClone(vafs=list(v1,v2), copyNumberCalls=list(cn1,cn2), sampleNames=c("Tumor1","tumor2"))
}
