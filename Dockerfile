FROM rocker/shiny:4.0.5

# System requirements
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    coinor-libcbc-dev coinor-libclp-dev libglpk-dev
    
## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean
    

ENV _R_SHLIB_STRIP_=true

RUN install2.r --error --skipinstalled \
shiny \
shinyjs \
dashboardthemes \
shinydashboard \
shinycssloaders \
shinyWidgets \
shinyBS \
shinyFiles \
highcharter \
ggvis \
png \
data.table \
DT \
optparse \
gprofiler2 \
heatmaply \
RSQLite \
DBI \
shinyjqui \
Seurat \
reticulate \
Matrix \
SparseM \
dplyr \
remotes \
BiocManager

RUN Rscript -e 'install.packages("plotly",version="4.9.4")'
RUN Rscript -e 'remotes::install_github("plotly/rasterly")'
RUN Rscript -e 'BiocManager::install("GenomicRanges")'

COPY /home/rstudio/scripts/SciViewerDev/pddExpn_0.1.0.tar.gz /srv/shiny-server/
COPY /home/rstudio/scripts/SciViewerDev/index.html /srv/shiny-server/
COPY /home/rstudio/scripts/SciViewerDev/README.md /srv/shiny-server/
COPY /home/rstudio/scripts/SciViewerDev/installs.R /srv/shiny-server/
COPY /home/rstudio/scripts/SciViewerDev/SciViewIn/* /srv/shiny-server/SciViewIn/
COPY /home/rstudio/scripts/SciViewerDev/SciView/* /srv/shiny-server/SciView/
COPY /home/rstudio/data/SCS/HS_healthyCholangiocytes_10x.h5ad /srv/shiny-server/data/

RUN R -e 'install.packages("/srv/shiny-server/pddExpn_0.1.0.tar.gz", repos = NULL, type = "source")'

USER shiny

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
