## server_compare.R

##-----------------------------------------------------------------
##    gene expression across available studies
##-----------------------------------------------------------------

observeEvent(c(compvariables$geneExpSummarized),{
  if(isTRUE(nrow(compvariables$geneExpSummarized)>0)){
    ## data table
    cat(" displaying summary expression table across all available studies for gene ",compvariables$gene,'\n')
    
    output$comparativeExpnTable <- reactable::renderReactable({
      myfun01 <- function(x){ log2(median(x,na.rm=T)+1) }
      qf <- compvariables$geneExp %>%
        dplyr::group_by(cell_type,feature,Database) %>%
        dplyr::summarize(MedExpr=myfun01(norm_avg_priorLT),expression_values = list(norm_avg_priorLT)) %>%
        mutate(boxplot = NA)
      qf$expression_values <- lapply(qf$expression_values, function(i){
        x <- ifelse(length(i)==1,paste0(rep(i,2),collapse = " "),paste0(i,collapse = " "))
        as.numeric(unlist(strsplit(x," ")))
      })
      reactable::reactable(qf, 
                defaultPageSize = 20,
                bordered = TRUE,
                defaultColDef = colDef(footerStyle = list(fontWeight = "italic",color='#41ab5d')),
                columns = list(
                  cell_type = colDef(filterable = T,searchable = T,sortable = T,resizable = T, footer = 'Cell type as annotated by author'),
                  feature = colDef(filterable = T,searchable = T,sortable = T,resizable = T,footer = "Gene name, as provided in the study"),
                  Database = colDef(filterable = T,searchable = T,sortable = T,resizable = T),
                  MedExpr = colDef(format = colFormat(digits = 2),footer = "Normalized median Expression Across donors",resizable = T),
                  expression_values = colDef(cell = function(values) {
                    sparkline(values, type = "bar", chartRangeMin = min(unlist(qf$expression_values)), chartRangeMax = max(unlist(qf$expression_values)))
                  },footer = 'Normalized expression per donor',resizable = T),
                  boxplot = colDef(cell = function(value, index) {
                    sparkline(qf$expression_values[[index]], type = "box")
                  },sortable = T,resizable = T,footer = 'Distribution of expression values')
                )
      )
    })
    
    #output$comparativeExpnTable <- DT::renderDataTable({
    #  pf <- compvariables$geneExpSummarized
    #  pf$AveExpr <- round(pf$AveExpr,2)
    #  pf <- pf[order(pf$AveExpr,decreasing = T),]
    #  datatable(pf,rownames = F,extensions = 'Buttons',
    #            options = list(searching = TRUE,pageLength = 10,dom = 'Bfrtip', 
    #                           buttons = list(list(extend = 'csv',  filename = paste0("ExpressionComparison_",compvariables$gene)),
    #                                          list(extend ='excel', filename = paste0("ExpressionComparison_",compvariables$gene))), 
    #                           scrollX = TRUE,scrollCollapse = TRUE)) %>% 
    #    formatStyle('AveExpr',
    #                background = styleColorBar(range(pf$AveExpr,na.rm=T), 'lightblue'),
    #                backgroundRepeat = 'no-repeat',backgroundPosition = 'left')
    #})
  }
})


observeEvent(c(input$comparativeExpnTable_rows_selected),{
  if(isTruthy(input$comparativeExpnTable_rows_selected) & isTRUE(nrow(compvariables$geneExpSummarized)>0)){
   
    output$comparativeExpnPlot <- renderUI({
      output$tmpl <- renderPlotly({
        db <- compvariables$geneExpSummarized[input$comparativeExpnTable_rows_selected,]$Database %>%
          as.character()
        pl <- compvariables$geneExp
        pl <- pl[pl$Database%in%db,]
        if(nrow(pl)>0){
          P <- pl %>%
            ggplot(aes(x=cell_type,y=norm_avg,fill=cell_type))+
            geom_boxplot(alpha=0.2) + facet_wrap(~Database,ncol=2,scales = 'free')+
            theme_bw()+
            theme(legend.position = 'none',strip.background = element_blank())+
            xlab("")+ylab("")+
            coord_flip()
          P %>% ggplotly()
        }else{
          empty_plot("select row from above table")
        }
      })
      plotlyOutput("tmpl")
    })
  }
})


##-----------------------------------------------------------------
##    gene summarized as marker across studies
##-----------------------------------------------------------------

observeEvent(c(compvariables$geneMarker),{
  if(isTRUE(nrow(compvariables$geneMarker)>0)){
    ## data table
    cat(" displaying summary of gene as marker across all available studies for gene ",compvariables$gene,'\n')
    output$comparativeMarkerTable <- DT::renderDataTable({
      cf <- compvariables$geneMarker
      cf$logFC = round(cf$logFC,2)
      cf$AveExpr = round(cf$AveExpr,2)
      cf$t = round(cf$t,2)
      cf$B = round(cf$B,2)
      cf <- cf[,c("Database","geneSymbol", "Test","logFC", "adj.P.Val")]
      datatable(cf,rownames = F,extensions = 'Buttons',
                options = list(searching = TRUE,pageLength = 5,dom = 'Bfrtip', 
                               buttons = list(list(extend = 'csv',  filename = paste0("Markers_FDR_CrossStudies")),
                                              list(extend ='excel', filename = paste0("Markers_FDR_CrossStudies"))), 
                               scrollX = TRUE,scrollCollapse = TRUE
                )) %>% 
        formatSignif(columns = c('adj.P.Val'),digits = 2) %>%
        formatStyle('logFC',
                    background = styleColorBar(range(cf$logFC,na.rm=T), 'lightblue'),
                    backgroundRepeat = 'no-repeat', backgroundPosition = 'right')
    })
  }else{
    output$comparativeMarkerTable <- DT::renderDataTable({c()})
  }
})







##-----------------------------------------------------------------
##    gene summarized as biological marker across studies
##-----------------------------------------------------------------


observeEvent(c(compvariables$geneBioMarker),{
  if(isTRUE(nrow(compvariables$geneBioMarker)>0)){
    ## data table
    cat(" displaying summary of gene as biomarker across all available studies for gene ",compvariables$gene,'\n')
    output$comparativeBioMarkerTable <- DT::renderDataTable({
      cf <- compvariables$geneBioMarker
      cf$logFC = round(cf$logFC,2)
      cf$AveExpr = round(cf$AveExpr,2)
      cf$t = round(cf$t,2)
      cf$B = round(cf$B,2)
      cf <- cf[,c("Database","geneSymbol", "Test","logFC", "adj.P.Val")]
      datatable(cf,rownames = F,extensions = 'Buttons',
                options = list(searching = TRUE,pageLength = 5,dom = 'Bfrtip', 
                               buttons = list(list(extend = 'csv',  filename = paste0("BioMarkers_FDR_CrossStudies")),
                                              list(extend ='excel', filename = paste0("BioMarkers_FDR_CrossStudies"))), 
                               scrollX = TRUE,scrollCollapse = TRUE
                )) %>% 
        formatSignif(columns = c('adj.P.Val'),digits = 2) %>%
        formatStyle('logFC',
                    background = styleColorBar(range(cf$logFC,na.rm=T), 'lightblue'),
                    backgroundRepeat = 'no-repeat', backgroundPosition = 'right')
    })
  }else{
    output$comparativeBioMarkerTable <- DT::renderDataTable({NULL})
  }
})

