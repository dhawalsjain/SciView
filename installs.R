
## CRAN
pckgs <- c("shiny","shinyjs","dashboardthemes","shinydashboard","shinycssloaders",
           "shinyWidgets","shinyBS","shinyFiles","highcharter","ggvis",
           "rasterly","png","data.table","DT","pddExpn","optparse","gprofiler2",
           "heatmaply","RSQLite","DBI","shinyjqui","Seurat","reticulate","Matrix",
           "SparseM","dplyr","BiocManager")
pckgs <- pckgs[!pckgs%in%installed.packages()]
install.packages(pckgs,dependencies = T,verbose = T)

## specific plotly version
install.packages("plotly", version='4.9.4')

## Bioconductor
pckgs <- c('GenomicRanges')
pckgs <- pckgs[!pckgs%in%installed.packages()]
BiocManager::install(pckgs)

## remote git
remotes::install_github("plotly/rasterly")

## local
install.packages("pddExpn_0.1.0.tar.gz", repos = NULL, type = "source")

