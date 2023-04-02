source("vars.R",local = TRUE,encoding = "UTF-8")
source("plot_functions.R",local = TRUE,encoding = "UTF-8")
source("MMdbLib.R",local = TRUE,encoding = "UTF-8")

##----------------------------------------------
warningMessage <- paste0("<font color=\"#e6550d\">
                           <br><br><br>
                           Requested plot/table is unavailable. <br>
                           Reasons:<br> 
                           1) Gene not present in the underlying database <br>
                           2) Gene locus is unavaialbe in Ensembl GRCh37 transcript database. <br>
                           3) A technical glitch with plot rendering? <br>
                           Try refreshing the page or contact the Data science team/ developers for details.
                           </font><br><br><br>")

##----------------------------------------------
js <- "
    function(el, x, inputName){
      var id = el.getAttribute('id');
      var gd = document.getElementById(id);
      var xScrollPos = el.scrollLeft || document.documentElement.scrollLeft;
      var yScrollPos = el.scrollTop || document.documentElement.scrollTop;
      var d3 = Plotly.d3;
      Plotly.plot(id).then(attach);
        function attach() {
          var xaxis = gd._fullLayout.xaxis;
          var yaxis = gd._fullLayout.yaxis;
          var l = gd._fullLayout.margin.l;
          var t = gd._fullLayout.margin.t;
          var coordinates = [null, null]
          
          gd.addEventListener('mousemove', 
          function(evt) {
              var coordinates = [xaxis.p2c(evt.pointerX-l), 
                                 yaxis.p2c(evt.pointerY-t)];
              Shiny.setInputValue(inputName, coordinates);
          });

        };
  }
  "

##----------------------------------------------

server <- function(input, output, session) {
  
  cat("DhawalJain\n")
  
  gc()
  Prog <- shiny::Progress$new()
  on.exit(Prog$close())
  Prog$set(message = "initializing the app..", value = 0)
  observe({ cat("tab item: ",input$AppTab,"\n") })

  ROOTS=c(workdir='.',
          datadir='C:/',
          home='/home/',
          shinydata='/srv/shiny-server/data/')
  
  ###---- 0) meta
  VARS <- reactiveValues(
    connGenes = list(),
    connList = list(),
    dbb = NULL,
    roots = ROOTS,
    sc_studies = c(),
    availGenes = NULL
  )
  
  variables <- reactiveValues(
    gene = NULL,
    availGenes=NULL,
    sc_study = NULL,
    sc_study_params = NULL,
    connID_1 = NULL,
    sc_study_attribs = NULL,
    donor = NULL,
    celltype = NULL,
    catvars = NULL,
    catvars_valuelist = NULL,
    contvars = NULL,
    availcelltypes = NULL,
    sc_louvain = NULL,
    sel_catvar = NULL,
    sel_catvarval = NULL,
    expndf=NULL,
    avail_markerTests = NULL,
    marker_cellTypeDf = NULL,
    marker_geneDf = NULL,
    cellmark_upenrich = list(),
    cellmark_dnenrich = list(),
    expndf_markers = NULL,
    avail_deTests = NULL,
    de_cellTypeDf = NULL,
    de_geneDf = NULL,
    cellde_upenrich = list(),
    cellde_dnenrich = list(),
    expndf_de = NULL
  ) 
  compvariables <- reactiveValues(
    gene = NULL,
    geneSymbol=NULL,
    geneExp = c(),
    geneMarker = c(),
    geneBioMarker = c(),
    geneExpSummarized = NULL
  )
  
  
  observeEvent({input$ab573dfg3},{
    if(isTruthy(input$ab573dfg3)){
      if(dir.exists(input$datadir)){
        VARS$roots = ROOTS
        inputdir=ifelse(isTRUE(grep("/$",input$datadir)),
                        as.character(input$datadir),
                        as.character(paste0(input$datadir,"/")))
        names(inputdir) = 'inputdir'
        cat("inputdir:",inputdir,"\n")
        VARS$roots = append(inputdir,VARS$roots)
        output$direxists <- NULL
        rm(inputdir)
      }else{
        cat('Directory does not exits!\n')
        output$direxists <- renderPrint({ cat('Directory does not exits!\n') })
      }
    }
  })
  observeEvent({VARS$roots},{
    if(isTruthy(VARS$roots)){
      shinyFileChoose(input,'infile', session=session,roots=VARS$roots, 
                      filetypes=c('', 'db'))
    }
  })
  
  ##--- make plot resizable and dragable
  #jqui_resizable(ui = "#sc_facetbyexpnviolin",operation = 'enable')
  #jqui_resizable(ui = "#sc_dedotplot",operation = 'enable' )
  #jqui_resizable(ui = "#sc_markerdotplot",operation = 'enable')
  #jqui_draggable(ui = "#sc_dedotplot,#sc_markerdotplot")
  
  ###---  1) genes database
  Prog$set(message = "gathering available studies..", value = 0.3)
  VARS$connGenes <- ({
    RSQLite::dbConnect(RSQLite::SQLite(),GENEFILE)
  })
  
  ###---- 2.0) ops on DB.txt file 
  output$availstudypaths <- renderText({
    xx <- tryCatch({
      xx <- read.delim(file = SCDBFILE,header = F,comment.char = "#",stringsAsFactors = F)$V1
      sapply(xx, function(u) file.exists(u))
    },error= function(e){
      NULL
    })
    paste0(xx,collapse = "\n")
  })
  observeEvent(c(input$ab573wer_dbgo,input$ab573wer_db),{
    if(isTruthy(input$ab573wer_dbgo) & isTruthy(input$ab573wer_db)){
      txt <- unlist(strsplit(input$ab573wer_db,","))
      txt <- unique(txt[txt != ""])
      txt <- sapply(txt, function(u) trimws(u, which = c("both"), whitespace = "[ \t\r\n]"))
      txt <- txt[sapply(txt, function(u) file.exists(u))]
      message(paste0("numbers of studies to append:",length(txt),"\n"))
      if(length(txt)>0){
        write.table(txt,file =SCDBFILE,append = T,quote = F,sep = "\n",row.names = F,col.names = F)
        session$reload()
      }
    }
  })
  observeEvent(c(input$ab573rty_dbgo,input$ab573rty_db),{
    if(isTruthy(input$ab573rty_dbgo) & isTruthy(input$ab573rty_db)){
      txt <- unlist(strsplit(input$ab573rty_db,","))
      txt <- unique(txt[txt != ""])
      txt <- sapply(txt, function(u) trimws(u, which = c("both"), whitespace = "[ \t\r\n]"))
      txt <- txt[sapply(txt, function(u) file.exists(u))]
      message(paste0("numbers of studies to replace:",length(txt),"\n"))
      if(length(txt)>0){
        write.table(txt,file =SCDBFILE,append = F,quote = F,sep = "\n",row.names = F,col.names = F)
        session$reload()
      }
    }
  })
  
  ###---- 2.1) set the sqlite databases with access 
  VARS$dbb <- ({
    tryCatch({
      xx <- read.delim(file = SCDBFILE,header = F,comment.char = "#",stringsAsFactors = F)$V1
      xx[sapply(xx, function(u) file.exists(u))]
    },error= function(e){
    })
  })
  observeEvent({input$infile},{
    VARS$dbb = NULL
    VARS$dbb = tryCatch({
      xx <- read.delim(file = SCDBFILE,header = F,comment.char = "#",stringsAsFactors = F)$V1
      xx[sapply(xx, function(u) file.exists(u))]
    },error= function(e){
    })
    cat(VARS$roots,"\n")
    inFile <- parseFilePaths(roots=VARS$roots, input$infile)
    if(isTruthy(inFile$datapath)){
      ext <- tools::file_ext(unname(inFile$datapath))
      validate(need(ext %in% c("db"), "Please select files with .db extension!"))
      VARS$dbb <- append(VARS$dbb,unname(inFile$datapath))
      VARS$dbb <- unique(VARS$dbb)
      rm(ext)
    }
  })
  
  ##--- 3) check if remote databases are allowed
  observe({
    if(USE_REMOTE_DB==TRUE){
      VARS$dbb = NULL
      cf <- GetStudyName(AppName=REPO_NAME)
      cf <- as.character(unlist(cf[,1]))
      cf <- cf[grep("^scRNA_",cf)]
      if(isTruthy(cf)){
        VARS$dbb <- unique(cf)
      }
      Prog$set(message = "Using remote databases..", value = 0.4)
      rm(cf)
    }else{
      Prog$set(message = "Using local databases..", value = 0.4)
      print("Currently working with local databases\n")
    }
  })
  
  
  ##--- 4) populate meta object
  observeEvent({VARS$dbb},{
    if(isTruthy(VARS$dbb)){
      VARS$connList <- ({
        connList <- list()
        if(length(VARS$dbb)>0){
          for(f in VARS$dbb){
            cat("populate meta object: ",f,"\n")
            id = pddExpn::randomStringGenerator(n = 1,lenght = 8)
            if(USE_REMOTE_DB==TRUE){
              connList[[id]] <- f #gsub("^scRNA_","",f)
            }else{
              connList[[id]] <- RSQLite::dbConnect(RSQLite::SQLite(),f)
            }
            rm(id,f)
          }
        }
        connList
      })
      VARS$sc_studies <- ({
        sc_studies <- c()
        for(i in 1:length(connList)){
          if(USE_REMOTE_DB==F){
            tablist <- RSQLite::dbListTables(connList[[i]])
            tablist <- tablist[grep("_study",tablist)]
          }else{
            tablist <- paste0(connList[[i]],"_study")
            #tablist <- gsub("^scRNA_","",tablist)
          }
          for(j in tablist){
            cat(i,"\t",j,"\n")
            query <- paste0("SELECT * FROM ",j)
            sc_studies <- rbind(sc_studies,
                                cbind(
                                  queryDB(HANDLER=connList[[i]], QUERY=query,
                                          REPO_NAME=REPO_NAME,
                                          USE_REMOTE_DB=USE_REMOTE_DB),
                                  ObjID=names(connList)[i]
                                )
            )
          }
          rm(i,j)
        }
        sc_studies
      })
    }
  })
  
  
  ##--- 5) list available studies
  if(T){
    output$browse_studytable <- DT::renderDataTable({
      if(isTruthy(nrow(VARS$sc_studies)>0)){
        cf <- VARS$sc_studies
        cf <- cf[,c("STATUS","ORGANISM","DISEASE","TISSUES","CellCount","SHORTDESCR")]
        names(cf) <- c("Status",'Species','Disease/Condition','Tissue','Number of Cells','Brief')
        cf$Species <- paste0(
          "<img src='",cf$Species,".png' height='30' data-toggle='tooltip' data-placement='right' title='",cf$Species,"'></img>"
        )
        cf$Status <- paste0(
          "<img src='",cf$Status,".png' height='30' data-toggle='tooltip' data-placement='right' title='",cf$Status,"'></img>"
        )
        datatable(cf,filter='top',autoHideNavigation = T, rownames = F,selection='single',extensions = 'Buttons',
                  options = list(pageLength = 20, autoWidth = TRUE),escape = FALSE) %>%
          formatStyle('Number of Cells',
                      background = styleColorBar(cf$'Number of Cells', '#a1d99b'),
                      backgroundSize = '95% 50%',
                      backgroundRepeat = 'no-repeat',
                      backgroundPosition = 'right')
      }
    })
  }
  
  if(F){
    output$browse_studytable <- reactable::renderReactable({
      if(isTruthy(nrow(VARS$sc_studies)>0)){
        cf <- VARS$sc_studies
        cf <- cf[,c("STATUS","ORGANISM","DISEASE","TISSUES","CellCount","SHORTDESCR")]
        names(cf) <- c("Status",'Species','Condition','Tissue','Number','Brief')
        reactable(cf,
                  defaultPageSize = 20,
                  bordered = TRUE,
                  defaultColDef = colDef(footerStyle = list(fontWeight = "italic",color='#41ab5d')),
                  columns = list(
                    Status = colDef(filterable = T,searchable = T,sortable = T,resizable = T, footer = 'Usage restrictions',cell = function(value) {
                      image <- img(src = sprintf("%s.png", value), style = "height: 24px;", alt = value)
                      tagList(div(style = "display: inline-block; width: 45px;", image),
                              value
                      )}),
                    Species = colDef(filterable = T,searchable = T,sortable = T,resizable = T, footer = 'Organism',cell = function(value) {
                      image <- img(src = sprintf("%s.png", value), style = "height: 24px;", alt = value)
                      tagList(div(style = "display: inline-block; width: 45px;", image),
                              value
                      )}),
                    Condition = colDef(name = "Biological Condition",filterable = T,searchable = T,sortable = T,resizable = T, footer = 'Disease/Experiment/Treatment-group etc'),
                    Tissue = colDef(filterable = T,searchable = T,sortable = T,resizable = T),
                    Number = colDef(name="Number of Cells",filterable = T,searchable = F,sortable = T,resizable = T,
                                    style = function(value) {
                                      bar_style(width = value / max(cf$Number), fill = "#2c5e77", color = "#fff")
                                    },
                                    align = "left",
                                    format = colFormat(digits = 1)),
                    Brief = colDef(name="Brief Study description",filterable = T,searchable = T,sortable = T,resizable = T)
                  )
        )
      }
    })
  }

  source(file = "server_df.R",local = TRUE,encoding = "UTF-8")$value
  Prog$set(message = "retrieved data for vizualisations..", value = 0.7)
  
  source(file = "server_sc.R",local = TRUE,encoding = "UTF-8")$value
  Prog$set(message = "prepared data for vizualisations..", value = 0.8)
  
  source(file = "server_compare.R",local = TRUE,encoding = "UTF-8")$value
  Prog$set(message = "retrieved data for cross study comparison..", value = 0.9)
  
  Prog$set(message = "done!", value = 1)
  
}