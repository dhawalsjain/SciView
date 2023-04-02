library(shiny)
library(shinyjs)
library(dashboardthemes)
library(shinydashboard)
library(shinycssloaders)
library(shinyWidgets)
library(shinyBS)
library(shinyFiles)
library(shinyjqui)

library(highcharter)
library(ggvis)
library(ggplot2)
library(ggridges)
library(plotly)
library(rasterly)
library(sparkline)
library(png)
library(data.table)
library(reactable)
library(sparkline)
library(DT)
library(pddExpn)
library(optparse)
library(gprofiler2)
library(heatmaply)

library(RSQLite)
library(DBI)
#library(RMySQL)

##-----------------------------------------------
GENEFILE="genes.db"
SCDBFILE="DBs.txt"
SOURCEREF='SciView'

##-- if using remote database
USE_REMOTE_DB = FALSE
REPO_NAME="Minuteman"

##--- Base expression parameter
BASE_EXPN_VAL=0


## Organism list for gprofiler2
ORGANISM <- list("Homo sapiens"='hsapiens',
                 "Homo sapien"='hsapiens',
                 "Rattus norvegicus" = "rnorvegicus",
                 "Mus musculus" = "mmusculus")




