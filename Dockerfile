# bioconductor 3.17
FROM bioconductor/bioconductor_docker@sha256:258c0e8b9d0e001adbca7f07372833e3f3b4c0495aa728502c6a2bb241383b18

RUN wget https://cran.r-project.org/src/contrib/Archive/NORMT3/NORMT3_1.0.4.tar.gz && \
    R CMD INSTALL NORMT3_1.0.4.tar.gz

RUN R --quiet -e "BiocManager::install(c('IRanges','MKmisc')); devtools::install_github('genome/bmm'); devtools::install_github('genome/sciClone')"

USER rstudio
