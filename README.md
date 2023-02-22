# SciView

Single Cell Interactive Viewer is an R shiny application that allows users to interactively visualize single cell datasets. The application contains two modules - 1) SciViewIn: The input application that processes single cell data from .h5ad or .seurat format into either a local or remote database. 2) SciView: The main application that allows interactive visualization. 

**Installation: From git**
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
- If running on ShinyServer, replace the **index.html** file with the provided copy. 

- Run the **SciViewIn** applications as following
```
runApp(appDir = "/path/to/SciViewIn/",launch.browser = T)
```
- Run the **SciView** applications as following
```
runApp(appDir = "/path/to/SciView/",launch.browser = T)
```

**Installation: From Docker**
- You can run SciView using docker image available through DockerHub 
- The contianer can be run either on and an AWS-EC2 instance or on local linux cluster (where docker is installed and working)
- You can listen to the SciView container on port 3838
- While running it on an AWS-EC2, please make sure that the port used for listening on DNS server is open for traffic and is unused by other applications.
- Make sure that docker is installed on the machine, before running following steps

- Pull Docker image
```
docker pull 24122016/sciview:latest
```
- Check port usage
```
lsof -i :8080 
lsof -i :3838 
```

- Remove other docker containers
```
docker rm -f $(docker ps -aq) 
```

- Doker docker container and and mount volume
      * In the below command, data under /home/rstudio/data directory of the EC2 instance is mounted into shiny server
      * The port 3838 from the docker container is exposed to 8080. Hence, one can listen the application on (DNS-server-address):8080
      * Mounted drive can be used to store the data (i.e. single cell seurat/h5ad ojects and .db files created by the application)

```
docker run -d -p 8080:3838 \
    -v /home:/srv/shiny-server/data \
    --name sciview sciview
```

- If you are using local linux machine for running SciView, you can listen to the application on 127.0.0.1:8080

- You can check the run logs on the above running docker instance by launching interactive shell

```
docker exec -it sciview bash

```








