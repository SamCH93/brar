## set R version (https://hub.docker.com/r/rocker/verse/tags)
FROM rocker/verse:4.5.0

## set up directories
RUN mkdir /home/rstudio/paper
COPY package /home/rstudio/package

## install R packages from CRAN the last day of the specified R version
RUN install2.r --error --skipinstalled --ncpus -1 \
    remotes knitr ggplot2 dplyr Ternary ggpubr ggh4x mvtnorm SimDesign && \
    R CMD INSTALL /home/rstudio/package/out/brar_0.1.tar.gz
    ## TODO install from GitHub with fixed commit hash
    ## R -e "remotes::install_github(repo = 'SamCH93/brar@ac4473dfcff753b1af2afa63168961b23ff4b2f1', subdir = 'package' upgrade_dependencies = FALSE)"
