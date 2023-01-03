library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinyWidgets)
library(shinyjqui)

library(Seurat)
library(reticulate)

library(data.table)
library(Matrix)
library(SparseM)
library(dplyr)

library(RSQLite)


ROOTS=c(workdir='.',
        datadir='C:/',
        home='/home/',
        shinydata='/srv/shiny-server/data/')

