tabItem(tabName = "scstudy",
        tags$hr(),
        fluidRow(column(width=3,
                        div(style='width:200px;overflow-x: auto;height:400px;overflow-y: auto;',
                            box(
                              selectInput("sel_celltype",label="celltypes",choices = NULL,multiple = T),
                              div(style="display: inline-block;vertical-align:top; horizontal-aling:center;",
                                  actionButton("sel_celltypego"," ",icon = icon("arrows-rotate"))
                              ),
                              solidHeader = T, collapsible = F,collapsed = F, width = 12
                            ),
                            box(
                              selectInput("sel_donor",label="donors",choices = NULL,multiple = T),
                              div(style="display: inline-block;vertical-align:top; horizontal-aling:center;",
                                  actionButton("sel_donorgo"," ",icon = icon("arrows-rotate"))
                              ),
                              solidHeader = T, collapsible = F,collapsed = F, width = 12
                            ),
                            box(
                              selectInput("sel_catvars",label="features",choices = NULL,multiple = F),
                              selectInput("sel_catvarsval",label="sub-features",choices = NULL,multiple = T),
                              div(style="display: inline-block;vertical-align:top; horizontal-aling:center;",
                                  actionButton("sel_metadatago"," ",icon = icon("arrows-rotate"))
                              ),
                              solidHeader = T, collapsible = F,collapsed = F, width = 12
                            )
                        ),
        ),
        column(width=6,
               plotlyOutput('umapplot'),
               #-- get plot 
               tags$br()
        ),
        column(width=3,
               fluidRow(column(width=12,
                               htmlOutput("umapplot_hover")),
                        style = "height:120px; background-color: white;"),
               tags$hr(),
               jqui_draggable(box(
                       selectizeInput("sel_gene",label="select gene",choices = NULL,multiple = T),
                       actionButton("sel_genego","Visualize!",icon = icon("arrows-rotate")),
                       solidHeader = T, collapsible = F,collapsed = F, width = 12
                       ),operation = "enable")
               
        )
        ),
        tags$hr(),
        tabBox(title = "",width = 14,side='left',
               tabPanel(title=" "),
               tabPanel(title="Data QC",
                        fluidRow(
                                column(width=12,
                                       verticalTabsetPanel(
                                               color = "gray",
                                               contentWidth = 11,
                                               verticalTabPanel(title=h6("Study description", style = 'font-size:13px;color:black;'),
                                                                icon=icon("info"),
                                                                box_height = "80px",
                                                                box_width = "80px",
                                                                htmlOutput("sc_study_attribs01")
                                               ),
                                               verticalTabPanel(title=h6("QC plots", style = 'font-size:13px;color:black;'),
                                                                box_height = "80px",
                                                                box_width = "80px",
                                                                icon=icon("eye-open", lib = "glyphicon"),
                                                                tags$style(HTML(".select-dropdown { font-size: 11px; }")),
                                                                fluidRow(column(width = 9,
                                                                                div(style='width:700px;overflow-x: auto;height:450px;overflow-y: auto;',
                                                                                    uiOutput("sc_qcscatter"))%>% shinycssloaders::withSpinner(color="#0dc5c1")
                                                                ),
                                                                column(width = 3,
                                                                       box(
                                                                               selectizeInput("qc_attrib02",label="Select QC attribute",choices = NULL,multiple = T),
                                                                               solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                                       )
                                                                )
                                                                ),
                                                                tags$br()
                                               )
                                       ),
                                       tags$br()
                                )
                        )
               ),
               tabPanel(title="Counts",
                        fluidRow(column(width = 12,
                                        verticalTabsetPanel(
                                                color = "gray",
                                                contentWidth = 11,
                                                verticalTabPanel(title=h6("Cell counts", style = 'font-size:13px;color:black;'),
                                                                 icon=icon("signal", lib = "glyphicon"),
                                                                 box_height = "100px",
                                                                 box_width = "80px",
                                                                 fluidRow(column(width = 9,
                                                                                 jqui_resizable(ui = highchartOutput("sc_qcbarplot") %>% shinycssloaders::withSpinner(color="#0dc5c1"),operation = 'enable')
                                                                                ),
                                                                          column(width = 3,
                                                                                box(
                                                                                 selectizeInput("qc_attrib01",label="Select attribute",choices = NULL,multiple = T),
                                                                                 radioButtons("rad_qcptgrp", label = shiny::HTML("<b>Aggregate by</b>"),
                                                                                             choices = list("Selected annotations" = 1, "Donors" = 2), 
                                                                                             selected = 1,inline = F),
                                                                                 conditionalPanel(
                                                                                         condition = "input.qc_attrib01.length>2",
                                                                                         selectizeInput("qc_attrib01a",label="Select sub-feature",choices = NULL,multiple = F)
                                                                                 ),
                                                                                 actionButton("sel_qcptgo","",icon = icon("arrows-rotate")),
                                                                                 solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                                                ),
                                                                                tags$br()
                                                                         )
                                                                 )
                                                )
                                   )
                                )
                        )
               ),
               tabPanel(title="Single-Gene",
                        fluidRow(
                                column(width=12,
                                       verticalTabsetPanel(
                                               color = "gray",
                                               contentWidth = 11,
                                               verticalTabPanel(title=h6("Scatter", style = 'font-size:13px;color:black;'),
                                                                icon=icon("signal", lib = "glyphicon"),
                                                                box_height = "80px",
                                                                box_width = "80px",
                                                                fluidRow(column(width=9, 
                                                                                #-- plot
                                                                                plotlyOutput('expnumapplot')
                                                                ),
                                                                column(width=3,
                                                                       box(
                                                                               selectInput("sc_gcolscale",label="select color-scale",choices = NULL,multiple = F),
                                                                               solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                                       )
                                                                )
                                                                ),
                                                                tags$br()
                                               ),
                                               verticalTabPanel(title=h6("Summary", style = 'font-size:13px;color:black;'),
                                                                icon=icon("signal", lib = "glyphicon"),
                                                                box_height = "80px",
                                                                box_width = "80px",
                                                                fluidRow(column(width = 9,
                                                                                highchartOutput("sc_gebarplot") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                                                                tags$br(),
                                                                                shiny::HTML("When expression is summarized per donor, only first attribute from the selected <b> Metadata feature</b> is used for summarizing the data.<br>")
                                                                ),
                                                                column(width = 3,
                                                                       box(
                                                                               selectInput("sc_expncatvar",label="Metadata feature",choices = NULL,multiple = T),
                                                                               radioButtons("rad_expnchoice", label = shiny::HTML("<b>Cells to use</b>"),
                                                                                            choices = list("All" = 1, "Only expressing" = 2), 
                                                                                            selected = 1,inline = F),
                                                                               radioButtons("rad_expngraphchoice", label = shiny::HTML("<b>Graph type</b>"),
                                                                                            choices = list("Barplot" = 1, "Boxplot (outlier omitted)" = 2), 
                                                                                            selected = 1,inline = F),
                                                                               radioButtons("rad_expnptgrp", label = shiny::HTML("<b>Group by</b>"),
                                                                                            choices = list("Cell" = 1, "Donor" = 2), 
                                                                                            selected = 1,inline = F),
                                                                               actionButton("sel_sgexpngo","",icon = icon("arrows-rotate")),
                                                                               solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                                       )
                                                                )
                                                                ),
                                                                tags$br()
                                               ),
                                               verticalTabPanel(title=h6("Counts", style = 'font-size:13px;color:black;'),
                                                                icon=icon("signal", lib = "glyphicon"),
                                                                box_height = "80px",
                                                                box_width = "80px",
                                                                fluidRow(column(width = 9,
                                                                                highchartOutput("sc_expncountbar") %>% shinycssloaders::withSpinner(color="#0dc5c1") ,
                                                                                tags$br(),
                                                                                shiny::HTML("The plot summarizes % of total cells expressing the genes along selected <b>Metadata feature</b><br>")
                                                                ),
                                                                column(width = 3,
                                                                       box(
                                                                               selectInput("sc_catvar_expncount",label="Metadata feature",choices = NULL,multiple = T),
                                                                               solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                                       )
                                                                )
                                                                ),
                                                                tags$br()
                                               )
                                       )
                                )
                        )
               ),
               tabPanel(title="Multi-Gene",
                        #shiny::HTML("Select multiple genes in the <b><i>Select gene</b></i> section above."),
                        fluidRow(column(width = 12,
                                        verticalTabsetPanel(
                                                color = "gray",
                                                contentWidth = 11,
                                                verticalTabPanel(title=h6("Features", style = 'font-size:13px;color:black;'),
                                                                 icon=icon("signal", lib = "glyphicon"),
                                                                 box_height = "80px",
                                                                 box_width = "80px",
                                                                 div(style='width:700px;overflow-x: auto;height:400px;overflow-y: auto;',
                                                                     uiOutput("sc_multigeneumap")),
                                                                 shiny::HTML("<b>NOTE:</b> The expression scale is arbritrary and the visualization should be used for qualitative comparison<br>"),
                                                                 tags$br()
                                                ),
                                                verticalTabPanel(title=h6("Dotplot", style = 'font-size:13px;color:black;'),
                                                                 icon=icon("signal", lib = "glyphicon"),
                                                                 box_height = "80px",
                                                                 box_width = "80px",
                                                                 fluidRow(column(width = 9,
                                                                                 #jqui_resizable(ui = uiOutput("sc_expndotplot"),operation = 'enable'),
                                                                                 div(style='width:700px;overflow-x: auto;height:400px;overflow-y: auto;',
                                                                                     uiOutput("sc_expndotplot")%>% shinycssloaders::withSpinner(color="#0dc5c1")),
                                                                                 shiny::HTML("The dotplot summarize average gene expression for selected genes.
                                                                                 Size of the point represent % of cells that express the gene. While color represents expression scale.<br>"), 
                                                                                 tags$br()
                                                                 ),
                                                                 column(width = 3,
                                                                        box(
                                                                                selectInput("sc_catvar_expndot",label="Metadata feature",choices = NULL,multiple = F),
                                                                                actionButton("sc_catvar_expndotgo","",icon = icon("arrows-rotate")),
                                                                                solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                                        )
                                                                 )
                                                                 )
                                                ),
                                                verticalTabPanel(title=h6("Correlations", style = 'font-size:13px;color:black;'),
                                                                 icon=icon("signal", lib = "glyphicon"),
                                                                 box_height = "80px",
                                                                 box_width = "80px",
                                                                 fluidRow(column(width = 6,
                                                                                 selectInput("sc_corrcelltype",label="select cell types",choices = NULL,multiple = T)
                                                                                 ),
                                                                          column(width = 3,
                                                                                 tags$br(),
                                                                                 actionButton("sc_corrcelltypego","",icon = icon("arrows-rotate"))
                                                                                 ),
                                                                          column(width = 3)
                                                                          ),
                                                                 shiny::HTML("<br>"),
                                                                 fluidRow(column(width = 6,
                                                                                 box(div(style='width:500px;overflow-x: auto;height:400px;overflow-y: auto;',
                                                                                         uiOutput("sc_corrplot") %>% shinycssloaders::withSpinner(color="#0dc5c1")),
                                                                                         solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                                                     )
                                                                                 ),
                                                                          column(width = 6,
                                                                                 box(div(style='width:500px;overflow-x: auto;height:400px;overflow-y: auto;',
                                                                                         uiOutput("sc_ridgeplot") %>% shinycssloaders::withSpinner(color="#0dc5c1")),
                                                                                     solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                                                     )
                                                                                 )
                                                                          ),
                                                                 box(fluidRow(column(width=3,
                                                                                     selectInput("sc_corrgenex",label="gene on x-axis",choices = NULL,multiple = F),
                                                                                     selectInput("sc_corrgeney",label="gene on y-axis",choices = NULL,multiple = F),
                                                                                     actionButton("sc_corrscattergo","",icon = icon("arrows-rotate"))
                                                                                     ),
                                                                              column(width=9,
                                                                                     jqui_resizable(uiOutput("sc_corrscatter"),operation = 'enable')
                                                                                     #uiOutput("sc_corrscatter") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                                                                                     )
                                                                              ),
                                                                     solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                                     ),
                                                                 shiny::HTML("<br>"),
                                                                 shiny::HTML("Pearson correlation coefficients are calculated between the genes along the cells expressing the genes.
                                                                     Based on the selection of cell types, the correlation coefficients are updated.<br>
                                                                     <br>The Ridgeplot summarizes expression distribution across the cells or selected group of cells.<br>
                                                                     <br><b>Note:</b> Minimum 2 genes are required for calculating correlation coefficients.<br><br><br>")
                                                ),
                                                verticalTabPanel(title=h6("FacetByExpn", style = 'font-size:13px;color:black;'),
                                                                 icon=icon("signal", lib = "glyphicon"),
                                                                 box_height = "80px",
                                                                 box_width = "80px",
                                                                 fluidRow(column(width = 3,
                                                                                 selectInput("sc_fbe",label="select base gene",choices = NULL,multiple = T)),
                                                                          column(width = 3,
                                                                                 tags$br(),
                                                                                 actionButton("sc_fbego","",icon = icon("arrows-rotate"))),
                                                                          column(width = 3),
                                                                          column(width = 3)
                                                                 ),
                                                                 tags$br(),
                                                                 fluidRow(column(width = 10,
                                                                                 jqui_resizable(ui = div(style='width:700px;overflow-x: hide;height:400px;overflow-y: hide;horizontal-aling:center;vertical-aling:center;',uiOutput("sc_facetbyexpnviolin") %>% shinycssloaders::withSpinner(color="#0dc5c1"),style="align:center"),
                                                                                                operation = 'enable')
                                                                 ),
                                                                 column(width = 2)
                                                                 ),
                                                                 tags$br(),
                                                                 shiny::HTML("The plot displays gene expressions after grouping the cells based on expression of the gene selection above. <br>
                                                                 Unlike correlation plots which are more quantitative and might miss subtle differences in expression patterns of the genes, this  
                                                                 visualization allows to qualitatively assess whether the pair of genes are correlated or not. <br><br><br>")
                                                ),
                                                verticalTabPanel(title=h6("CoExpn", style = 'font-size:13px;color:black;'),
                                                                 icon=icon("signal", lib = "glyphicon"),
                                                                 box_height = "80px",
                                                                 box_width = "80px",
                                                                 fluidRow(column(width = 3),
                                                                          column(width = 6,
                                                                                 div(style='width:300px;overflow-x: hide;height:300px;overflow-y: hide;',
                                                                                     uiOutput("sc_coexpnumap") %>% shinycssloaders::withSpinner(color="#0dc5c1"))
                                                                          ),
                                                                          column(width = 3)
                                                                         ),
                                                                 tags$br(),
                                                                 fluidRow(column(width = 2),
                                                                          column(width = 8,
                                                                                 div(style='width:700px;overflow-x: hide;height:400px;overflow-y: hide;horizontal-aling:center;vertical-aling:center;',
                                                                                     uiOutput("sc_coexpnbar") %>% shinycssloaders::withSpinner(color="#0dc5c1"),style="align:center")
                                                                          ),
                                                                          column(width = 2)
                                                                          ),
                                                                 tags$br(),
                                                                 shiny::HTML("Gene expression is defined as non-zero expression of the gene. Genes selected in the above section are used to display co-expression.<br>
                                                                                 Cells that express <i><u>all genes</i></u> are displayed in <font color=\"#3f007d\"><b> dark-blue </font></b> color, while those that either express some or none are colored as <font color=\"#efedf5\"><b> light-blue </font></b>.<br>
                                                                                 The cell-counts are further summarized in adjacent barplot.<br>
                                                                                             <b>Note:</b> Minimum 2 genes are required for displaying co-expression.<br><br><br>")
                                                )
                                                
                                        )
                        )
                )
               ),
               tabPanel(title="Cell Markers",
                        fluidRow(column(width = 12,
                                        verticalTabsetPanel(
                                                color = "gray",
                                                contentWidth = 11,
                                                verticalTabPanel(title=h6("In selected genes", style = 'font-size:13px;color:black;'),
                                                                 icon=icon("signal", lib = "glyphicon"),
                                                                 #box_height = "80px",
                                                                 #box_width = "80px",
                                                                 shiny::HTML("<u>Select genes from the above section and <b> hit the refresh button.</b></u><br><br>"),
                                                                 jqui_resizable(ui = div(style='width:700px;overflow-x: auto;height:450px;overflow-y: auto;',
                                                                                         uiOutput("sc_markerdotplot")),
                                                                                operation = 'enable'),
                                                                 shiny::HTML("<i>blank plot suggests that selected genes do not qualify as marker at FDR 20%</i>")
                                                ),
                                                verticalTabPanel(title=h6("In specific comparisons", style = 'font-size:13px;color:black;'),
                                                                 icon=icon("signal", lib = "glyphicon"),
                                                                 #box_height = "80px",
                                                                 #box_width = "80px",
                                                                 fluidRow(column(width = 6,
                                                                                 selectInput("sel_markercelltype",label="select celltype",choices = NULL,multiple = F)
                                                                                 ),
                                                                          column(width = 6,
                                                                                 selectInput("sel_markercomp",label="select comparison",choices = NULL,multiple = F)
                                                                                 ),
                                                                          style = "height:80px; background-color: white;"
                                                                 ),
                                                                 fluidRow(column(width = 6,
                                                                                 sliderInput("sel_markerfdrslider", label = shiny::HTML("<b>FDR cutoff</b>"), min=0,max=0.2, value = 0.05,step = 0.05)
                                                                                 ),
                                                                          column(width = 6,
                                                                                 box(solidHeader = T, collapsible = F,collapsed = F, width = 12,
                                                                                     actionButton("sc_markeropgo",shiny::HTML("calculate enrichments!"),icon = icon("arrows-rotate")),
                                                                                     shiny::HTML("<br><i>Enrichment results displayed in the next section</i>")
                                                                                     )
                                                                                 ),
                                                                          style = "height:80px; background-color: white;"
                                                                 ),
                                                                 shiny::HTML("<br><br>"),
                                                                 fluidRow(column(width=6,
                                                                                 shiny::HTML("<b>Marker table</b><br>Select one or more rows from the table below to generate gene expression plots.<br>"),
                                                                                 DT::dataTableOutput("marker_tab") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                                                                                 ),
                                                                          column(width=6,
                                                                                 shiny::HTML("<b>Volcano plot</b><br>The volcano plot summarizes top 500 genes sorted by their FDR values<br>"),
                                                                                 plotlyOutput("sc_markerVolcano") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                                                                                 )#,
                                                                          #style = "height:250px; background-color: white;"
                                                                 ),
                                                                 shiny::HTML("<br><br>"),
                                                                 div(style='width:700px;overflow-x: auto;height:450px;overflow-y: auto;', uiOutput("sc_markerumap"))
                                                                 #shiny::HTML("<br><br>Enrichments for the the comparison selected here are displayed in the next section. Hit <u>calculate enrichment</u><br>")
                                                ),
                                                verticalTabPanel(title=h6("Enrichments", style = 'font-size:13px;color:black;'),
                                                                 icon=icon("signal", lib = "glyphicon"),
                                                                 fluidRow(column(width = 6,
                                                                                 radioButtons("sel_markerdirchange", label = shiny::HTML("<b>Direction of change</b>"),
                                                                                              choices = list("Over-expressed" = 1, "Under-expressed" = 2), 
                                                                                              selected = 1,inline = T)
                                                                                 ),
                                                                          column(width=6,
                                                                                 shiny::HTML("This view is linked to the FDR selection made on previous section - <b> Specific comparison</b>. Enrichments and the displays are produced using <a href=\"https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html\",target=\"_blank\"> gProfiler2 </a>.</i><br><br>")
                                                                                 )
                                                                          ),
                                                                 shiny::HTML("<br><br>"),
                                                                 uiOutput('markerEnrichmentText'),
                                                                 shiny::HTML("<br><br>"),
                                                                 uiOutput("enrichmarkers_plot")%>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                                                 DT::dataTableOutput("enrichmarker_tab"),
                                                                 shiny::HTML("<br><br>")
                                                )
                                        )
                                        )
                        )
               ),
               tabPanel(title="Biological Markers",
                        fluidRow(column(width=12,
                                        verticalTabsetPanel(
                                                color = "gray",
                                                contentWidth = 11,
                                                verticalTabPanel(title=h6("In\nSelected\ngenes", style = 'font-size:13px;color:black;'),
                                                                 icon=icon("signal", lib = "glyphicon"),
                                                                 #box_height = "120px",
                                                                 #box_width = "80px",
                                                                 fluidRow(column(width = 12,
                                                                                 conditionalPanel(
                                                                                         condition = "input.sel_genego>'0'",
                                                                                         selectInput("sel_detags",label="select celltype",choices = NULL,multiple = F)
                                                                                 )
                                                                                 )
                                                                          ),
                                                                 #shiny::HTML("<u>Select genes from the above section and <b> hit the refresh button.</b></u><br><br>"),
                                                                 shiny::HTML("<br>"),
                                                                 jqui_resizable(ui = div(style='width:700px;overflow-x: auto;height:450px;overflow-y: auto;',
                                                                                         uiOutput("sc_dedotplot")),
                                                                                operation = 'enable'),
                                                                 shiny::HTML("<i>blank plot suggests that selected genes do not qualify as biological markers at FDR 20%</i>")
                                                ),
                                                verticalTabPanel(title=h6("In\nSpecific\ncomparisons", style = 'font-size:13px;color:black;'),
                                                                 icon=icon("signal", lib = "glyphicon"),
                                                                 #box_height = "120px",
                                                                 #box_width = "80px",
                                                                 fluidRow(column(width = 6,
                                                                                 selectInput("sel_decelltype",label="select celltype",choices = NULL,multiple = F)
                                                                                 ),
                                                                          column(width = 6,
                                                                                 selectInput("sel_decomp",label="select comparison",choices = NULL,multiple = F)
                                                                                 ),
                                                                          style = "height:80px; background-color: white;"
                                                                          ),
                                                                 fluidRow(column(width = 6,
                                                                                 sliderInput("sel_defdrslider", label = shiny::HTML("<b>FDR cutoff</b>"), min=0,max=0.2, value = 0.05,step = 0.05)
                                                                                 ),
                                                                          column(width = 6,
                                                                                 box(solidHeader = T, collapsible = F,collapsed = F, width = 12,
                                                                                     actionButton("sc_deopgo",shiny::HTML("calculate enrichments!"),icon = icon("arrows-rotate")),
                                                                                     shiny::HTML("<br><i>Enrichment results displayed in the next section</i>")
                                                                                     )
                                                                                 ),
                                                                          style = "height:80px; background-color: white;"
                                                                          ),
                                                                 shiny::HTML("<br><br>"),
                                                                 fluidRow(column(width=6,
                                                                                 shiny::HTML("<b>Biological Marker table</b><br>Select one or more rows from the table below to generate gene expression plots.<br>"),
                                                                                 DT::dataTableOutput("de_tab") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                                                                                 ),
                                                                          column(width=6,
                                                                                 shiny::HTML("<b>Volcano plot</b><br>The volcano plot summarizes top 500 genes sorted by their FDR values<br>"),
                                                                                 plotlyOutput("sc_deVolcano") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                                                                                 )#,
                                                                          #style = "height:250px; background-color: white;"
                                                                          ),
                                                                 shiny::HTML("<br><br>"),
                                                                 div(style='width:700px;overflow-x: auto;height:450px;overflow-y: auto;',uiOutput("sc_deumap"))
                                                                 #shiny::HTML("<br><br>Enrichments for the the comparison selected here are displayed in the next section. Hit <u>calculate enrichment</u><br>")
                                                ),
                                                verticalTabPanel(title=h6("Enrichments", style = 'font-size:13px;color:black;'),
                                                                 icon=icon("signal", lib = "glyphicon"),
                                                                 fluidRow(column(width = 6,
                                                                                 radioButtons("sel_dedirchange", label = shiny::HTML("<b>Direction of change</b>"),
                                                                                              choices = list("Over-expressed" = 1, "Under-expressed" = 2), 
                                                                                              selected = 1,inline = T)
                                                                                 ),
                                                                          column(width=6,
                                                                                 shiny::HTML("This view is linked to the FDR selection made on previous section - <b> Specific comparison</b>. Enrichments and the displays are produced using <a href=\"https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html\",target=\"_blank\"> gProfiler2 </a>.</i><br><br>")
                                                                                 )
                                                                          ),
                                                                 shiny::HTML("<br><br>"),
                                                                 uiOutput('deEnrichmentText'),
                                                                 shiny::HTML("<br><br>"),
                                                                 uiOutput("enrichde_plot")%>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                                                 shiny::HTML("<br><br>"),
                                                                 DT::dataTableOutput("enrichde_tab")
                                                )
                                        )
                                )
                        )

                )## add here
        )
)
