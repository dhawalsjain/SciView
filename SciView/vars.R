library(shiny)
library(shinyjs)
library(dashboardthemes)
library(shinydashboard)
library(shinycssloaders)
library(shinyWidgets)
library(shinyBS)
library(shinyFiles)

library(highcharter)
library(ggvis)
library(plotly)
library(rasterly)
library(png)
library(data.table)
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

## Organism list for gprofiler2
ORGANISM <- list("Homo sapiens"='hsapiens',
                 "Rattus norvegicus" = "rnorvegicus",
                 "Mus musculus" = "mmusculus")




