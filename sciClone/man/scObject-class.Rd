\name{scObject-class}
\Rdversion{1.0}
\docType{class}
\alias{scObject-class}

\title{Class "scObject" ~~~ }

\description{
Objects of this class hold data used by the sciClone package, 
including variant vaf data, copy number calls, clustering results,
density and purity calculations, etc.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("scObject")}.
}

\section{Slots}{
  \describe{
    \item{\code{vafs.merged}:}{Object of class \code{"data.frame"}
       a data frame containing variant positions, vafs, cn, clusters, etc
    }
    \item{\code{clust}:}{Object of class \code{"list"}
      a list containing the results of clustering the vafs
    }
    \item{\code{densities}:}{Object of class \code{"list"}
      list containing information on kernel density and peaks
    }
    \item{\code{dimensions}:}{Object of class \code{"numeric"}
      the number of dimensions that are being used for clustering
    }
    \item{\code{marginalClust}:}{Object of class \code{"list"}
      marginal clustering data
    }
    \item{\code{sampleNames}:}{Object of class \code{"character"}
      a vector of names corresponding to the input samples
    }
    \item{\code{vafs.1d}:}{Object of class \code{"list"}
      a list containing 1d vafs and clusters for each sample
    }
    \item{\code{ellipse.metadata}:}{Object of class \code{"list"}
      metadata related to confidence metric plotting
    }
  }
}
\section{Methods}{
No methods defined with class "scObject" in the signature.
}
\keyword{classes}
