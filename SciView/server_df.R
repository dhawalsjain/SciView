## server_df.R

##-----------------------------------------------
##--- 1) study details, available genes etc
##-----------------------------------------------
observeEvent({input$browse_studytable_rows_selected},{
  if(isTruthy(input$browse_studytable_rows_selected)){
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = "resetting/initializing dataframe extractions..", value = 0)
    cat(" Initializing the data\n")
    
    variables$gene = NULL
    variables$availGenes=NULL
    variables$sc_study = NULL
    variables$sc_study_params = NULL
    variables$connID_1 = NULL
    variables$sc_study_attribs = NULL
    variables$donor = NULL
    variables$celltype = NULL
    variables$catvars = NULL
    variables$catvars_valuelist = NULL
    variables$contvars = NULL
    variables$availcelltypes = NULL
    variables$sc_louvain = NULL
    variables$sel_catvar = NULL
    variables$sel_catvarval = NULL
    variables$expndf = NULL
    variables$avail_markerTests = NULL
    variables$marker_cellTypeDf = NULL
    variables$marker_geneDf = NULL
    variables$cellmark_upenrich = list()
    variables$cellmark_dnenrich = list()
    variables$expndf_markers = NULL
    variables$avail_deTests = NULL
    variables$de_cellTypeDf = NULL
    variables$de_geneDf = NULL
    variables$cellde_upenrich = list()
    variables$cellde_dnenrich = list()
    variables$expndf_de = NULL

    compvariables$gene = input$sel_genex
    compvariables$geneExp = c()
    compvariables$geneMarker = c()
    compvariables$geneBioMarker = c()
    compvariables$geneExpSummarized = NULL
    
    #-- get the study data frame
    cf <- VARS$sc_studies[input$browse_studytable_rows_selected,]
    variables$sc_study_params = cf
    
    #-- description
    variables$sc_study_attribs <- ({
      paste("<font color=\"#FF0000\"><b> STUDY STATUS: </b> ",as.character(cf$STATUS),"! </font><br>",
            "<b> Species: </b>", as.character(cf$ORGANISM),"<br>",
            "<b> Number of Cells: </b>", format(as.numeric(cf$CellCount),big.mark=",")," ",
            "<b> Number of features: </b>", format(as.numeric(cf$FeatureCount),big.mark=","),"<br>",
            "<b> TISSUE: </b>", as.character(cf$TISSUES)," ",
            "<b> DISEASE/CONDITION: </b>", as.character(cf$DISEASE)," ",
            "<b> SAMPLE SIZE: </b>", as.character(cf$SampleSize),"<br>",
            "<b> PUBMED: </b>", a(as.character(cf$PMID),href=paste0("https://pubmed.ncbi.nlm.nih.gov/",as.character(cf$PMID)), target="_blank"),
            "<b> ACCSSION: </b>", a("Data",href=paste0(as.character(cf$GEO)), target="_blank"),"<br>",
            "<b> STUDY ABSTRACT: </b> <font color=\"#bdbdbd\">", as.character(cf$Description),"</font><br>"
      )
    })
    
    #--- study
    variables$sc_study = as.character(cf$Database)
    
    #--- connID
    variables$connID_1 <- ({as.character(cf$ObjID) })
    
    #--- base dataframe
    variables$sc_louvain <- ({
      query <- paste0("SELECT * FROM ",variables$sc_study,"_metaFeatures")
      louvain <- queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
                         QUERY=query,REPO_NAME=REPO_NAME,
                         USE_REMOTE_DB=USE_REMOTE_DB)
      louvain$value <- NA
      louvain
    })
    scProg$set(message = "retrieved required base dataframe..", value = 0.5)
    
    #--- celltype
    x <- unlist(strsplit(as.character(cf$cell_type),","))
    if(!is.null(x)){
      variables$celltype=x
    }else{
      scProg$set(message = "can't locate user-defined cell type column in the database. Exiting!", value = 0.5)
      stop("can't locate user-defined cell type column in the database. Exiting!")
    }
    rm(x)
    
    #--- donor_var
    x <- unlist(strsplit(as.character(cf$donor),","))
    if(!is.null(x)){
      variables$donor=x
    }
    
    #--- categorical variables
    x <- unlist(strsplit(as.character(cf$CATVAR),","))
    if(!is.null(variables$donor)){x <- x[x!=variables$donor]}
    x <- x[x!=variables$celltype]
    if(length(x)>1){
      variables$catvars=x
    }
    rm(x)
    
    #-- call available unique values for catvars
    if(!is.null(variables$catvars)){
      for(f in variables$catvars){
        if(f %in% names(variables$sc_louvain)){
          variables$catvars_valuelist[[f]] <- unique(variables$sc_louvain[,f])
        }else{
          variables$catvars <- variables$catvars[variables$catvars!=f]
        }
      }
    }

    #--- continuous variables
    x <- unlist(strsplit(as.character(cf$CONTVAR),","))
    if(!is.null(x)){
      variables$contvars=x
    }
    rm(x)
    
    #-- available cell types
    variables$availcelltypes <- ({
      unique(as.character(variables$sc_louvain[,as.character(variables$celltype)]))
    })
    
    #-- genes
    cat(" getting gene lists for: ", as.character(cf$Database), " objID: ", as.character(cf$ObjID),"\n")
    query <- paste0("SELECT geneSymbol FROM ",as.character(cf$Database),"_data")
    availGenes <- queryDB(HANDLER=VARS$connList[[as.character(cf$ObjID)]], 
                          QUERY=query,REPO_NAME=REPO_NAME,
                          USE_REMOTE_DB=USE_REMOTE_DB)
    availGenes <- tryCatch({
      query <- paste0("SELECT * FROM genehomologs WHERE species = '",as.character(cf$ORGANISM),"'")
      g1 <- RSQLite::dbGetQuery(VARS$connGenes, query)
      g1 <- g1[,c("geneSymbol","DBClassKey")] %>% unique()
      availGenes <- merge(availGenes,g1,by="geneSymbol",all.x=T)
      g1 <- g1[-grep("^ENS",g1$geneSymbol),]
      names(g1)[1] <- "Name"
      availGenes <- merge(availGenes,g1,by="DBClassKey",all.x=T)
      availGenes$Name <- ifelse(is.na(availGenes$Name),as.character(availGenes$geneSymbol),as.character(availGenes$Name))
      availGenes <- availGenes[,c(2,3,1)]
      availGenes
    },error=function(e){
      availGenes <- cbind(availGenes,Name=availGenes$geneSymbol,DBClassKey=NA)
      availGenes
    })
    updateSelectizeInput(session, 'sel_gene', choices = unique(availGenes[,2]), 
                         server = TRUE,selected ="ITGB8")
    updateSelectizeInput(session, 'sel_genex', choices = unique(availGenes[,2]), 
                         server = TRUE,selected ="ITGB8")
    variables$availGenes = availGenes
    
    #--- gene expression dataframe initialization
    variables$expndf = data.frame(SAMPID=variables$sc_louvain[,'SAMPID'])
    rownames(variables$expndf) = variables$expndf[,1]
    variables$expndf[,1] = NULL
    variables$expndf_markers = variables$expndf
    variables$expndf_de = variables$expndf
    
    scProg$set(message = "done..", value = 1)
  }
})

##-----------------------------------------------
##--- 2) alias gene names
##-----------------------------------------------
if(F){
  gAlias <- reactive({
    query <- paste0("SELECT * FROM gAliasProtein")
    gGenes <- RSQLite::dbGetQuery(VARS$connGenes, query)
    gAlias <- list()
    for(f in input$sel_gene){
      hgnc <- gGenes[gGenes$geneSymbol==f,]$HGNC %>% as.character
      gg <- as.character(gGenes[gGenes$HGNC%in%hgnc,]$geneSymbol)
      gAlias[[f]] <- unique(c(f,gg))
    }
    gAlias
  })
}

if(T){
  gAlias <- reactive({
    gAlias <- list()
    gg <- variables$availGenes
    for(f in input$sel_gene){
      gAlias[[f]] <- unique(c(f,gg[gg$Name==f,]$geneSymbol[1]))
      cat("Gene Alias:",f,": ",paste0(gAlias[[f]],collapse = ","),"\n")
    }
    gAlias
  })

}

observeEvent({variables$sc_study_attribs},{
  if(isTruthy(variables$sc_study_attribs)){
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = "retrieving alias for gene names, if any..", value = 0)
    #-- messages
    scProg$set(message = paste0(" study1: ", variables$sc_study), value = 0.1)
    cat(" study1: ", variables$sc_study,"\n")
    scProg$set(message = paste0(" study1 connID: ", variables$connID_1), value = 0.2)
    cat(" study1 connID: ", variables$connID_1,"\n")
    scProg$set(message = paste0(" cell type column: ",variables$celltype), value = 0.3)
    cat(" cell type column: ",variables$celltype,"\n")
    scProg$set(message = paste0(" Donor column: ",variables$donor), value = 0.4)
    cat(" Donor column: ",variables$donor,"\n")
    scProg$set(message = paste0(" Categorical variables: ",variables$catvars), value = 0.5)
    cat(" Categorical variables: ",variables$catvars,"\n")
    scProg$set(message = paste0(" Continuous variables: ",variables$contvars), value = 0.6)
    cat(" Continuous variables: ",variables$contvars,"\n")
    scProg$set(message = paste0(" Available celltypes:  ",variables$availcelltypes), value = 0.7)
    cat(" Available celltypes:  ",variables$availcelltypes,"\n")
    scProg$set(message = paste0(" Alias gene names: ", paste0(unlist(gAlias()),collapse = ";")), value = 0.8)
    cat(" Alias gene names: ", paste0(unlist(gAlias()),collapse = ";"),"\n")
    cat("\n")
    scProg$set(message = "done..", value = 1)
  }
})


##-----------------------------------------------
##--- 3) update selection choices
##-----------------------------------------------



##-----------------------------------------------
##--- 4) gene expression dataframes
##-----------------------------------------------
observeEvent({input$sel_genego},{
    if(isTruthy(variables$sc_louvain) & isTruthy(input$sel_genego) &  isTruthy(input$sel_gene)){
      if(is.null(variables$gene) | sum(variables$gene!=input$sel_gene)>0 ){
        scProg <- shiny::Progress$new()
        on.exit(scProg$close())
        scProg$set(message = "retrieving plot data for multiple genes..", value = 0)
        cat(" Generating dataframe list for multiple genes..\n")
        
        #gg <- variables$availGenes
        #variables$gene=gg[gg$Name%in%input$sel_gene,]$geneSymbol  
        #cat("variable genes: ",paste0(variables$gene),"\n")
        variables$gene=input$sel_gene  
        
        ##------------ generate and update expndf
        for(g in input$sel_gene){
          y <- g
          if(!y %in% names(variables$expndf)){ #names(variables$pldf)
            pl <- NULL
            cntr =1 
            while(cntr<=length(gAlias()[[g]]) && !isTRUE(nrow(pl)>0)){
              cat("  gene/alias name: ", gAlias()[[g]][cntr],"\n")
              pl <-  get_plot_df_sc04(VARS$connList[[variables$connID_1]],
                                      variables$sc_study,
                                      gAlias()[[g]][cntr],
                                      variables$sc_louvain,
                                      REPO_NAME=REPO_NAME,
                                      USE_REMOTE_DB=USE_REMOTE_DB)
              cntr = cntr +1 
            }
            if(isTRUE(nrow(pl)>0)){
              scProg$set(message = paste0('  getting plot df for gene ',g, ":",y,":",names(variables$expndf)), value = 0.2)
              cat('  getting plot df for gene ',g, ":",y,":",names(variables$expndf),"\n")
              variables$expndf <- cbind(variables$expndf,pl$value)
              names(variables$expndf)[ncol(variables$expndf)] <- g
            }
          }else{
            scProg$set(message = paste0('plot dataframe for gene ',g, " is already calculated"), value = 0.2)
            cat(' plot dataframe for gene ',g, " is already calculated\n")
          }
        }
        rm(y,cntr,g)
        z <- paste0(input$sel_gene)
        for(y in names(variables$expndf)){ 
          if(!y%in%z){
            scProg$set(message = paste0('  removing plot dataframe for gene ',y), value = 0.5)
            cat('  removing plot dataframe for gene ',y, "\n")
            if(y %in%names(variables$expndf)){
              variables$expndf[,y] <- NULL
            }
          }
        }
        scProg$set(message = "done..", value = 1)
      }
    }
  })


##-----------------------------------------------
##--- 5) marker gene dataframes
##-----------------------------------------------
observeEvent(c(variables$gene),{
  if(isTruthy(variables$gene) & isTruthy(input$sel_genego)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("querying selected gene as marker across available features"), value = 0)
    cat("querying selected gene as marker across available features\n")
    variables$marker_geneDf = NULL
    
    cf <- local({
      cf <- c()
      for(g in variables$gene){
        pl <- c()
        cntr=1
        scProg$set(message = paste0("querying gene for celltype marker: ",g), value = 0.5)
        cat("querying gene for celltype marker: ",g,"\n")
        while(cntr<=length(gAlias()[[g]]) && !isTRUE(nrow(pl)>0)){
          pl <- tryCatch({
            query <- paste0("SELECT * FROM ",variables$sc_study,"_Marker WHERE geneSymbol = '",gAlias()[[g]][cntr],"'")
            queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
                    QUERY=query,REPO_NAME=REPO_NAME,
                    USE_REMOTE_DB=USE_REMOTE_DB)
          },error=function(e){
            NULL
          })
          cntr = cntr + 1
        }
        if(isTRUE(nrow(pl)>0)){
          cf <- rbind(cf,pl)
        }
        rm(pl)
      }
      cf
    })
    if(isTRUE(nrow(cf)>0)){
      scProg$set(message = paste0("marker dataframe for ", variables$gene," has ",unique(cf$Test)," unique tests"), value = 0.7)
      cat("marker dataframe for ", variables$gene," has ",unique(cf$Test)," unique tests","\n")
      variables$marker_geneDf <- cf
    }
    rm(cf)
    gc()
    scProg$set(message = "done..", value = 1)
  }
})

observeEvent(c(input$sel_markercelltype),{
  if(isTruthy(input$sel_markercelltype) & isTruthy(variables$sc_study)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("retrieving marker genes for ",input$sel_markercelltype), value = 0)
    cat("retrieving marker genes for ",input$sel_markercelltype,"\n")
    #--- resetting the dataframes
    variables$avail_markerTests = NULL
    variables$marker_cellTypeDf = NULL
    variables$expndf_markers = NULL
    variables$expndf_markers = data.frame(SAMPID=variables$sc_louvain[,'SAMPID'])
    rownames(variables$expndf_markers) = variables$expndf_markers[,1]
    variables$expndf_markers[,1] = NULL

    query <- paste0("SELECT * FROM ",variables$sc_study,"_Marker WHERE Tag = '",input$sel_markercelltype,"'")
    pl <- tryCatch({
      queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
              QUERY=query,REPO_NAME=REPO_NAME,
              USE_REMOTE_DB=USE_REMOTE_DB)
    },error=function(e){
      return(NULL)
    })
    if(isTRUE(nrow(pl)>0)){
      cat("  marker dataframe for ", input$sel_markercelltype,": ",nrow(pl),"genes; test: ",unique(pl$Test),"\n")
      variables$avail_markerTests = unique(pl$Test)
      variables$marker_cellTypeDf <- pl
    }
    rm(pl)
    gc()
    scProg$set(message = "done..", value = 1)
  }
})

observeEvent(c(input$marker_tab_rows_selected,input$sel_markerfdrslider,input$sel_markercomp),{
  if(isTruthy(input$marker_tab_rows_selected) & isTruthy(variables$marker_cellTypeDf) &
     isTruthy(variables$sc_study) & isTruthy(variables$sc_louvain) &
     isTruthy(input$sel_markercomp) & isTruthy(input$sel_markerfdrslider)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("extracting expression values for selected marker genes from the table"), value = 0)
    cat("extracting expression values for selected marker genes from the table\n")
    cat("  ..row index: ",input$marker_tab_rows_selected,"\n")
    
    cf <- variables$marker_cellTypeDf
    cf <- cf[as.numeric(cf$adj.P.Val) <= as.numeric(input$sel_markerfdrslider),]
    cf <- cf[cf$Test==input$sel_markercomp,]
    cf <- cf[order(cf$adj.P.Val),]
    
    ID = cf[input$marker_tab_rows_selected,]$geneSymbol %>%
      as.character()
    scProg$set(message = paste0("marker gene/s ",ID), value = 0.2)
    cat("marker gene/s ",ID,"...\n")
    ##------------ generate and update expression data frame
    for(g in ID){
      if(!g %in% names(variables$expndf_markers)){
        pl <- NULL
        pl <- get_plot_df_sc04(VARS$connList[[variables$connID_1]],
                               variables$sc_study,
                               g,
                               variables$sc_louvain,
                               REPO_NAME=REPO_NAME,
                               USE_REMOTE_DB=USE_REMOTE_DB)
        if(isTRUE(nrow(pl)>0)){
          cat('  getting plot df for gene ',g,"\n")
          variables$expndf_markers <- cbind(variables$expndf_markers,pl$value)
          names(variables$expndf_markers)[ncol(variables$expndf_markers)] <- g
        }
      }else{
        scProg$set(message = paste0("expression entry for marker gene ",g, " is already extracted"), value = 0.5)
        cat('expression entry for marker gene ',g, " is already extracted\n")
      }
    }
    rm(y,g)
    z <- ID
    for(y in names(variables$expndf_markers)){ 
      if(!y%in%z){
        scProg$set(message = paste0("removing expression entry for selected marker gene ",y), value = 0.5)
        cat('removing expression entry for selected marker gene ',y, "\n")
        if(y %in%names(variables$expndf_markers)){
          variables$expndf_markers[,y] <- NULL
        }
      }
    }
    scProg$set(message = "done..", value = 1)
  }
})

observeEvent(c(input$sc_markeropgo),{
  if(isTruthy(input$sc_markeropgo) & isTruthy(variables$sc_study) & 
     isTruthy(input$sel_markercelltype) &  isTRUE(nrow(variables$marker_cellTypeDf)>0) & 
     isTruthy(input$sel_markerfdrslider) & isTruthy(input$sel_markercomp)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("computing enrichment of the marker genes"), value = 0)
    cat("computing enrichment of the marker genes\n")
    cat("Organism: ", variables$sc_study_params[,"ORGANISM"])

    
    variables$cellmark_upenrich = list()
    variables$cellmark_dnenrich = list()
    pl <- variables$marker_cellTypeDf 
    pl <- pl[pl$Test==input$sel_markercomp,]
    ## Enrichment
    if(isTRUE(nrow(pl)>0)){
      variables$cellmark_upenrich <- local({
        qf <-  pl
        qf <- qf[as.numeric(qf$adj.P.Val) <= as.numeric(input$sel_markerfdrslider) & qf$logFC>0,]$geneSymbol
        qf <- unique(qf)
        if(length(qf)>2){
          qf = list(A=qf)
          names(qf) <- 'Upregulated (wrt rest cells)'
          XX = gprofiler2::gost(query = qf[1],organism =  ORGANISM[[variables$sc_study_params[,"ORGANISM"]]],
                           ordered_query = F,exclude_iea = T,
                           significant = F,user_threshold = 0.05,correction_method = 'g_SCS')
          XX[['info']]  <-     shiny::HTML("<u><b>Celltype: </b>",input$sel_markercelltype,"       
                 <b>Comparison: </b>",input$sel_markercomp,"   
                 <b>FDR: </b>",input$sel_markerfdrslider,"            
                 <b>Direction of change: </b>Over-expressed </u><br>")
          cat(" XX names:" ,paste0(names(XX),collapse = ","),"\n")
          XX
        }else{ NULL }
      })
      scProg$set(message = paste0("calculated enrichment for upregulated genes"), value = 0.5)
      cat("calculated enrichment for upregulated genes\n")
      variables$cellmark_dnenrich <- local({
        qf <-  pl
        qf <- qf[as.numeric(qf$adj.P.Val) <= as.numeric(input$sel_markerfdrslider) & qf$logFC < 0,]$geneSymbol
        qf <- unique(qf)
        if(length(qf)>2){
          qf = list(A=qf)
          names(qf) <- 'Downregulated (wrt rest cells)'
          XX <- gprofiler2::gost(query = qf[1],organism = ORGANISM[[variables$sc_study_params[,"ORGANISM"]]],#ORGANISM[[variables$sc_study[,"ORGANISM"]]],#
                           ordered_query = F,exclude_iea = T,
                           significant = F,user_threshold = 0.05,correction_method = 'g_SCS')
          XX[['info']]  <-     shiny::HTML("<u><b>Celltype: </b>",input$sel_markercelltype,"       
                 <b>Comparison: </b>",input$sel_markercomp,"  
                 <b>FDR: </b>",input$sel_markerfdrslider,"            
                 <b>Direction of change: </b>Under-expressed </u><br>")
          XX
        }else{ NULL }
      })
      scProg$set(message = paste0("calculated enrichment for downregulated genes"), value = 0.9)
      cat("calculated enrichment for downregulated genes\n")
      gc()
      scProg$set(message = "done..", value = 1)
    }
  }
})



##-----------------------------------------------
##--- 6) biological marker gene dataframes
##-----------------------------------------------
observeEvent(c(variables$gene),{
  if(isTruthy(variables$gene) & isTruthy(input$sel_genego)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("querying selected gene as biological marker across available comparisons"), value = 0)
    cat("querying selected gene as biological marker across available comparisons\n")
    variables$de_geneDf = NULL
    
    cf <- local({
      cf <- c()
      for(g in variables$gene){
        pl <- c()
        cntr=1
        scProg$set(message = paste0("querying gene for biological marker: ",g), value = 0.5)
        cat("querying gene for bioloical marker: ",g,"\n")
        while(cntr<=length(gAlias()[[g]]) && !isTRUE(nrow(pl)>0)){
          pl <- tryCatch({
            query <- paste0("SELECT * FROM ",variables$sc_study,"_DEG WHERE geneSymbol = '",gAlias()[[g]][cntr],"'")
            queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
                    QUERY=query,REPO_NAME=REPO_NAME,
                    USE_REMOTE_DB=USE_REMOTE_DB)
          },error=function(e){
            NULL
          })
          cntr = cntr+1
        }
        if(isTRUE(nrow(pl)>0)){
          cf <- rbind(cf,pl)
        }
        rm(pl)
      }
      cf
    })
    
    if(isTRUE(nrow(cf)>0)){
      scProg$set(message = paste0("biological marker dataframe for ", variables$gene," has ",unique(cf$Test)," unique tests"), value = 0.7)
      cat("biological marker dataframe for ", variables$gene," has ",unique(cf$Test)," unique tests","\n")
      #cf$Test <- paste0(cf$Test,"\n(",cf$geneSymbol,")")
      variables$de_geneDf <- cf
    }
    rm(cf)
    gc()
    scProg$set(message = "done..", value = 1)
  }
})

observeEvent(c(input$sel_decelltype),{
  if(isTruthy(input$sel_decelltype) & isTruthy(variables$sc_study)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("retrieving biological marker genes for ",input$sel_decelltype), value = 0)
    cat("retrieving biological marker genes for ",input$sel_decelltype,"\n")
    #--- resetting the dataframes
    variables$avail_deTests = NULL
    variables$de_cellTypeDf = NULL
    variables$expndf_de = NULL
    variables$expndf_de = data.frame(SAMPID=variables$sc_louvain[,'SAMPID'])
    rownames(variables$expndf_de) = variables$expndf_de[,1]
    variables$expndf_de[,1] = NULL
    
    query <- paste0("SELECT * FROM ",variables$sc_study,"_DEG WHERE Tag = '",input$sel_decelltype,"'")
    pl <- tryCatch({
      queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
              QUERY=query,REPO_NAME=REPO_NAME,
              USE_REMOTE_DB=USE_REMOTE_DB)
    },error=function(e){
      return(NULL)
    })
    if(isTRUE(nrow(pl)>0)){
      cat("Biological marker dataframe for ", input$sel_decelltype,": ",nrow(pl),"genes; test: ",unique(pl$Test),"\n")
      variables$avail_deTests = unique(pl$Test)
      variables$de_cellTypeDf <- pl
    }
    rm(pl)
    gc()
    scProg$set(message = "done..", value = 1)
  }
})

observeEvent(c(input$de_tab_rows_selected),{
  if(isTruthy(input$de_tab_rows_selected) & isTruthy(variables$de_cellTypeDf) &
     isTruthy(variables$sc_study) & isTruthy(variables$sc_louvain)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("extracting expression values for selected biological marker genes from the table"), value = 0)
    cat("extracting expression values for selected biological marker genes from the table\n")
    cat(" .. row index: ",input$de_tab_rows_selected,"\n")
    
    cf <- variables$de_cellTypeDf
    cf <- cf[as.numeric(cf$adj.P.Val) <= as.numeric(input$sel_defdrslider),]
    cf <- cf[cf$Test==input$sel_decomp,]
    cf <- cf[order(cf$adj.P.Val),]
    
    ID = cf[input$de_tab_rows_selected,]$geneSymbol %>%
      as.character()
    scProg$set(message = paste0("marker gene/s ",ID), value = 0.2)
    cat("marker gene/s ",ID,"...\n")
    ##------------ generate and update expression data frame
    for(g in ID){
      if(!g %in% names(variables$expndf_de)){
        pl <- NULL
        pl <- get_plot_df_sc04(VARS$connList[[variables$connID_1]],
                               variables$sc_study,
                               g,
                               variables$sc_louvain,
                               REPO_NAME=REPO_NAME,
                               USE_REMOTE_DB=USE_REMOTE_DB)
        if(isTRUE(nrow(pl)>0)){
          cat('  getting plot df for gene ',g,"\n")
          variables$expndf_de <- cbind(variables$expndf_de,pl$value)
          names(variables$expndf_de)[ncol(variables$expndf_de)] <- g
        }
      }else{
        scProg$set(message = paste0("expression entry for biological marker gene ",g, " is already extracted"), value = 0.5)
        cat('expression entry for biological marker gene ',g, " is already extracted\n")
      }
    }
    rm(y,g)
    z <- ID
    for(y in names(variables$expndf_de)){ 
      if(!y%in%z){
        scProg$set(message = paste0("removing expression entry for selected biological marker gene ",y), value = 0.5)
        cat('removing expression entry for selected biological marker gene ',y, "\n")
        if(y %in%names(variables$expndf_de)){
          variables$expndf_de[,y] <- NULL
        }
      }
    }
    scProg$set(message = "done..", value = 1)
  }
})

observeEvent(c(input$sc_deopgo),{
  if(isTruthy(input$sc_deopgo) & isTruthy(variables$sc_study_params) & 
     isTruthy(input$sel_decelltype) &  isTRUE(nrow(variables$de_cellTypeDf)>0) & 
     isTruthy(input$sel_defdrslider) & isTruthy(input$sel_decomp)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("computing enrichment of the biological marker genes"), value = 0)
    cat("computing enrichment of the biological marker genes\n")
    
    pl <- variables$de_cellTypeDf 
    pl <- pl[pl$Test==input$sel_decomp,]
    
    variables$cellde_upenrich = list()
    variables$cellde_dnenrich = list()
    
    ## Enrichment
    if(isTRUE(nrow(pl)>0)){
      variables$cellde_upenrich <- local({
        qf <-  pl
        qf <- qf[as.numeric(qf$adj.P.Val) <= as.numeric(input$sel_defdrslider) & qf$logFC>0,]$geneSymbol
        qf <- unique(qf)
        if(length(qf)>2){
          qf = list(A=qf)
          names(qf) <- 'Upregulated (wrt rest cells)'
          XX = gprofiler2::gost(query = qf[1],organism = ORGANISM[[variables$sc_study_params[,"ORGANISM"]]],#'hsapiens',#
                                ordered_query = F,exclude_iea = T,
                                significant = F,user_threshold = 0.05,correction_method = 'g_SCS')
          XX[['info']]  <-     shiny::HTML("<u><b>Celltype: </b>",input$sel_decelltype,"       
                 <b>Comparison: </b>",input$sel_decomp,"   
                 <b>FDR: </b>",input$sel_defdrslider,"            
                 <b>Direction of change: </b>Over-expressed </u><br>")
          XX
        }else{ NULL }
      })
      scProg$set(message = paste0("calculated enrichment for upregulated genes"), value = 0.5)
      cat("calculated enrichment for upregulated genes\n")
      
      variables$cellde_dnenrich <- local({
        qf <-  pl
        qf <- qf[as.numeric(qf$adj.P.Val) <= as.numeric(input$sel_defdrslider) & qf$logFC < 0,]$geneSymbol
        qf <- unique(qf)
        #cat("... organism:",ORGANISM[[variables$sc_study_params[,"ORGANISM"]]],"\n")
        if(length(qf)>2){
          qf = list(A=qf)
          names(qf) <- 'Downregulated (wrt rest cells)'
          XX <- gprofiler2::gost(query = qf[1],organism = ORGANISM[[variables$sc_study_params[,"ORGANISM"]]], #'hsapiens',#
                                 ordered_query = F,exclude_iea = T,
                                 significant = F,user_threshold = 0.05,correction_method = 'g_SCS')
          XX[['info']]  <-     shiny::HTML("<u><b>Celltype: </b>",input$sel_decelltype,"       
                 <b>Comparison: </b>",input$sel_decomp,"  
                 <b>FDR: </b>",input$sel_defdrslider,"            
                 <b>Direction of change: </b>Under-expressed </u><br>")
          XX
        }else{ list() }
      })
      scProg$set(message = paste0("calculated enrichment for downregulated genes"), value = 0.9)
      cat("calculated enrichment for downregulated genes\n")
      
      gc()
      scProg$set(message = "done..", value = 1)
    }
  }
})


##-----------------------------------------------
##--- 7) gene across studies
##-----------------------------------------------
observeEvent(c(input$sel_genexgo),{
  if(input$AppTab=='compare' & isTruthy(input$sel_genexgo) & 
     isTruthy(input$sel_genex) & isTruthy(VARS$sc_studies)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("extracting expression of ",input$sel_genex,"across studies"), value = 0)
    message(paste0("extracting expression of ",input$sel_genex,"across studies\n"))

    compvariables$gene = input$sel_genex
    compvariables$geneExp = c()
    compvariables$geneMarker = c()
    compvariables$geneBioMarker = c()
    compvariables$geneExpSummarized = NULL
    
    gAliasX <- local({
      gAliasX <- list()
      gg <- variables$availGenes
      dbclasskey = unique(gg[gg$Name==compvariables$gene,]$DBClassKey)[1] ## only use the 1st 
      query <- paste0("SELECT * FROM genehomologs WHERE DBClassKey = '",dbclasskey,"'")
      gg <- RSQLite::dbGetQuery(VARS$connGenes, query)
      gg <- gg[,c("geneSymbol","species")] %>% unique()
      gAliasX <- split(gg$geneSymbol,f = gg$species)
      gAliasX
    })
    scProg$set(message = paste0(" .. retrieved gene aliases"), value = 0.1)
    message(paste0(" ..retrieved gene aliases\n"))
    
    for(i in 1:nrow(VARS$sc_studies)){
      pl <- NULL
      ql <- NULL
      gl <- NULL
      cntr =1 
      cf <- VARS$sc_studies[i,]
      cat(" extracting data for the gene",compvariables$gene," in ",cf$Database,"\n")
      gAliasXs <- gAliasX[cf$ORGANISM]
      if(length(gAliasXs)>0){
        while(cntr<=length(gAliasXs[[1]]) && !isTRUE(nrow(pl)>0)){ ##1 because we are only querying one gene; no multiple genes 
          cat("  gene/alias name: ", gAliasXs[[1]][cntr],"\n")
          ##-- 01) Expression
          scProg$set(message = paste0(" .. retrieving expression from ",cf$Database), value = i/nrow(VARS$sc_studies))
          message(paste0(" .. retrieving expression from ",cf$Database,"\n"))
          pl <- local({
            tryCatch({
              query <- paste0("SELECT * FROM ",as.character(cf$Database),"_FeatureSummaryKeys")
              x <- queryDB(HANDLER=VARS$connList[[as.character(cf$ObjID)]], 
                      QUERY=query,REPO_NAME=REPO_NAME,
                      USE_REMOTE_DB=USE_REMOTE_DB)
              query <- paste0("SELECT * FROM ",as.character(cf$Database),"_FeatureSummary WHERE feature = '",gAliasXs[[1]][cntr],"'")
              y <- queryDB(HANDLER=VARS$connList[[as.character(cf$ObjID)]], 
                      QUERY=query,REPO_NAME=REPO_NAME,
                      USE_REMOTE_DB=USE_REMOTE_DB)
              x$norm_avg_priorLT <- unname(unlist(y[,2:ncol(y)]))
              x$feature <- unname(unlist(unique(y[1,1])))
              message(paste0(names(x),collapse = ","))
              x
            },error=function(e){ NULL
            })
          })
          ##-- 02) Marker
          scProg$set(message = paste0(" .. retrieving marker from ",cf$Database), value = i/nrow(VARS$sc_studies))
          message(paste0(" .. retrieving marker from ",cf$Database,"\n"))
          gl <- local({
            query <- paste0("SELECT * FROM ",as.character(cf$Database),"_Marker WHERE geneSymbol = '",gAliasXs[[1]][cntr],"'")
            tryCatch({
              queryDB(HANDLER=VARS$connList[[as.character(cf$ObjID)]], 
                      QUERY=query,REPO_NAME=REPO_NAME,
                      USE_REMOTE_DB=USE_REMOTE_DB)
            },error=function(e){ NULL
            })
          })
          ##-- 03) BioMarker
          scProg$set(message = paste0(" .. retrieving disease marker from ",cf$Database), value = i/nrow(VARS$sc_studies))
          message(paste0(" .. retrieving disease marker from ",cf$Database,"\n"))
          ql <- local({
            query <- paste0("SELECT * FROM ",as.character(cf$Database),"_DEG WHERE geneSymbol = '",gAliasXs[[1]][cntr],"'")
            tryCatch({
              queryDB(HANDLER=VARS$connList[[as.character(cf$ObjID)]], 
                      QUERY=query,REPO_NAME=REPO_NAME,
                      USE_REMOTE_DB=USE_REMOTE_DB)
            },error=function(e){ NULL
            })
          })
          ## counter
          cntr = cntr +1 
        }
      }
      if(isTRUE(nrow(pl)>0)){
        pl$Database <- as.character(cf$Database)
        compvariables$geneExp <- rbind(compvariables$geneExp, pl)
        cat("Expression Comparison: A study added\n")
      }
      if(isTRUE(nrow(gl)>0)){
        gl$Database <- as.character(cf$Database)
        compvariables$geneMarker <- rbind(compvariables$geneMarker, gl)
        cat("Marker Expression Comparison: A study added\n")
      }
      if(isTRUE(nrow(ql)>0)){
        ql$Database <- as.character(cf$Database)
        compvariables$geneBioMarker <- rbind(compvariables$geneBioMarker, ql)
        cat("Biological Marker Expression Comparison: A study added\n")
      }
      rm(pl,ql,gl)
      gc()
    }
    
    ## summarized table
    myfun01 <- function(x){ log2(mean(x,na.rm=T)+1) }
    compvariables$geneExpSummarized <- ({
      compvariables$geneExp %>%
        dplyr::group_by(cell_type,feature,Database) %>%
        dplyr::summarize(AveExpr=myfun01(norm_avg_priorLT))
    })
    scProg$set(message = paste0(" .. done!"), value = 1)
    

    
  }
})






