An R package for inferring the subclonal architecture of tumors

## Installation instructions:

As of mid-2022, the NORMT3 package, which is a dependency of sciclone/bmm, has been removed from CRAN. It can be installed manually by doing something like: 

    $ wget https://cran.r-project.org/src/contrib/Archive/NORMT3/NORMT3_1.0.4.tar.gz
    $ R CMD install NORMT3_1.0.4.tar.gz
    
Then proceed with the below instructions:

I forked sciclone from the orginial repo to fix the 'Error in xtfrm.data.frame(x) : cannot xtfrm data frames' issue in the sciClone() function. Make sure dplyr is installed (e.g. install.packages('tidyverse') and R >= 4.1.0.

Both the 'sciClone' package and it's 'bmm' dependency can be installed by doing the following:

    #install IRanges from bioconductor
    source("http://bioconductor.org/biocLite.R")
    biocLite("IRanges")
    #install devtools if you don't have it already
    install.packages("devtools")
    library(devtools)
    install_github("genome/bmm")
    #install_github("genome/sciClone")
    # install sciClone with fix for 'Error in xtfrm.data.frame(x) : cannot xtfrm data frames'
    install_github("kunstner/sciClone")

If you prefer to build the package by hand, follow these steps:

- Make sure that you have the dependencies from the CRAN and BioConductor repos:
IRanges, rgl, RColorBrewer, ggplot2, grid, plotrix, methods, NORMT3, MKmisc, TeachingDemos, dplyr

- install the bmm package from [https://github.com/genome/bmm](https://github.com/genome/bmm)

- Download and build from source:

        # git clone git@github.com:genome/sciclone.git
        git clone git@github.com:kunstner/sciclone.git
        R CMD build sciclone
        R CMD INSTALL sciClone_1.1.0.tar.gz

## Usage
    library(sciClone)

    #read in vaf data from three related tumors
    #format is 5 column, tab delimited: 
    #chr, pos, ref_reads, var_reads, vaf

    v1 = read.table("data/vafs.tumor1.dat",header=T);
    v2 = read.table("data/vafs.tumor2.dat",header=T);
    v3 = read.table("data/vafs.tumor3.dat",header=T);

    #read in regions to exclude (commonly LOH)
    #format is 3-col bed
    regions = read.table("data/exclude.loh")

    #read in segmented copy number data
    #4 columns - chr, start, stop, segment_mean   
    cn1 = read.table("data/copy_number_tum1")
    cn2 = read.table("data/copy_number_tum2")
    cn3 = read.table("data/copy_number_tum3")

    #set sample names
    names = c("Sample1","Sample2","Sample3")


    #Examples:
    #------------------------------------
    #1d clustering on just one sample
    sc = sciClone(vafs=v1,
             copyNumberCalls=cn1,
             sampleNames=names[1],
             regionsToExclude=reg1)
    #create output
    writeClusterTable(sc, "results/clusters1")
    sc.plot1d(sc,"results/clusters1.1d.pdf")

    #------------------------------------
    #2d clustering using two samples:
    sc = sciClone(vafs=list(v1,v2),
                  copyNumberCalls=list(cn1,cn2),
                  sampleNames=names[1:2],
                   regionsToExclude=regions)
    #create output
    writeClusterTable(sc, "results/clusters2")
    sc.plot1d(sc,"results/clusters2.1d.pdf")
    sc.plot2d(sc,"results/clusters2.2d.pdf")


    #------------------------------------
    #3d clustering using three samples:
    sc = sciClone(vafs=list(v1,v2,v3),
                  copyNumberCalls=list(cn1,cn2,cn3),
                  sampleNames=names[1:3],
                   regionsToExclude=regions)
    #create output
    writeClusterTable(sc, "results/clusters2")
    sc.plot1d(sc,"results/clusters2.1d.pdf")
    sc.plot2d(sc,"results/clusters2.2d.pdf")
    sc.plot3d(sc, sc@sampleNames, size=700, outputFile="results/clusters3.3d.gif")

    #This pattern generalizes up to N samples, except for plotting, which caps out at 3d for obvious reasons.

## Visualization

### single-tumor plot
![1d plot](http://i.imgur.com/n4JNs9t.png)

### 2d comparison plot
![2d plot](http://i.imgur.com/8h0qAWx.png)

### 3d comparison plot
![3d plot](http://i.imgur.com/iM0V1kq.gif)

## Notes

- Requires host system to have imagemagick installed before it can produce animated gif output of 3d plots.

- Input formats described in more detail in the R documentation (see `?sciClone`)

- Many questions regarding sciClone usage have been asked and answered on Biostar: https://www.biostars.org/t/sciclone/

## Accessory Scripts and Data
The [sciClone-meta](https://github.com/genome/sciclone-meta) repo contains all data and scripts used to create the figures in the manuscript. It also contains a small suite of tests that demonstrate the capabilities of sciClone and verify that it is installed correctly.

## Reference
Manuscript published at [PLoS Computational Biology (doi:10.1371/journal.pcbi.1003665)](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003665)

SciClone: Inferring clonal architecture and tracking the spatial and temporal patterns of tumor evolution

Christopher A. Miller<sup>1\*</sup>, Brian S. White<sup>2\*</sup>, Nathan D. Dees<sup>1</sup>, John S. Welch<sup>2,3</sup>, Malachi Griffith<sup>1</sup>, Obi Griffith<sup>1</sup>, Ravi Vij<sup>2,3</sup>, Michael H. Tomasson<sup>2,3</sup>, Timothy A. Graubert<sup>2,3</sup>, Matthew J. Walter<sup>2,3</sup>, William Schierding<sup>1</sup>, Timothy J. Ley<sup>1,2,3</sup>, John F. DiPersio<sup>2,3</sup>, Elaine R. Mardis<sup>1,3,4</sup>, Richard K. Wilson<sup>1,3,4</sup>, and Li Ding<sup>1,2,3,4</sup>

<sup>1</sup>The Genome Institute

<sup>2</sup>Department of Medicine

<sup>3</sup>Siteman Cancer Center

<sup>4</sup>Department of Genetics Washington University, St. Louis, MO 63110, USA

<sup>*</sup> These authors contributed equally to this work
