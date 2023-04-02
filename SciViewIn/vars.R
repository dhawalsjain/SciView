library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinyWidgets)
library(shinyjqui)

library(Seurat)
library(reticulate)
library(anndata)

library(data.table)
library(Matrix)
library(SparseM)
library(dplyr)

library(RSQLite)

library(future)
library(parallel)

ROOTS=c(workdir='.',
        datadir='C:/',
        home='/home/',
        shinydata='/srv/shiny-server/data/',
        dhawal_local='/home/rstudio/data/SCS/')

MCORES=round(parallel::detectCores(all.tests = FALSE, logical = TRUE)/3)
options(future.globals.maxSize= 1048576000)
