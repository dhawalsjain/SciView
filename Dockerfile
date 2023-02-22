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
    coinor-libcbc-dev coinor-libclp-dev libglpk-dev \
    r-base \
    python3.8 \
    python3-dev \
    python3-pip \
    python3-venv

CMD alias python=/usr/bin/python3.8
CMD alias python3=/usr/bin/python3.8
RUN ln -s /usr/bin/python3.8 /usr/bin/python
    
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
BiocManager \
reshape \
ggthemes

RUN Rscript -e 'install.packages("plotly",version="4.9.4")'
RUN Rscript -e 'remotes::install_github("plotly/rasterly")'
RUN Rscript -e 'BiocManager::install("GenomicRanges")'
RUN Rscript -e 'install.packages("anndata",version="0.7.5.3")'

## Data and Scripts
COPY pddExpn_0.1.0.tar.gz /srv/shiny-server/ 
COPY index.html /srv/shiny-server/ 
COPY README.md /srv/shiny-server/ 
COPY installs.R /srv/shiny-server/ 
COPY Dockerfile /srv/shiny-server/ 
COPY SciViewIn/ /srv/shiny-server/SciViewIn/ 
COPY SciView/ /srv/shiny-server/SciView/ 
RUN mkdir -p /srv/shiny-server/data

RUN Rscript -e 'install.packages("/srv/shiny-server/pddExpn_0.1.0.tar.gz", repos = NULL, type = "source")'

USER shiny

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
