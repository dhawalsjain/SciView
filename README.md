# SciView

Single Cell Interactive Viewer is an R shiny application that allows users to interactively visualize single cell datasets. The application contains two modules - 1) SciViewIn: The input application that processes single cell data from .h5ad or .seurat format into either a local or remote database. 2) SciView: The main application that allows interactive visualization. 

**Installation**
- The application is developed in R 4.0 and uses several dependancies listed below. 
- R Shiny libraries
    * library(shiny)
    * library(shinydashboard)
    * library(shinyFiles)
    * library(shinyWidgets)
    * library(shinyjs)
    * library(dashboardthemes)
    * library(shinycssloaders)
    * library(shinyBS)
- R data handling libraries
    * library(Seurat)
    * library(reticulate)
    * library(data.table)
    * library(Matrix)
    * library(SparseM)
    * library(dplyr)
    * library(preprocessCore)
    * library(DT)
    * library(optparse)
    * library(gprofiler2)
- R data visualization libraries
    * library(highcharter)
    * library(plotly) #version='4.9.4'
    * library(rasterly)
    * library(png)
    * library(heatmaply)
- R database libraries
    * library(RSQLite)
    * library(RMySQL)
    * library(DBI)
- Custom R library
    * library(pddExpn) : Available in this repo

- Set up R locally or on your R-server. 
- Clone the repo and install all above-listed dependencies.
- Copy **genes.db** data file from https://drive.google.com/file/d/1GFiNSQK5KiZ7U0ytv6Lbkod_ln4ZiifD/view?usp=share_link and store under the folder **SciView**

- Run the **SciViewIn** applications as following
```
runApp(appDir = "/path/to/SciViewIn/",launch.browser = T)
```
- Run the **SciView** applications as following
```
runApp(appDir = "/path/to/SciView/",launch.browser = T)
```


