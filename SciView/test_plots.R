####------------------------------------------------------------------------------
####----------- test dataframe limits
####------------------------------------------------------------------------------
if(F){
  pl <- c()
  
  cf <- data.frame(x=1,y=1)
  
  
  cf$y = paste0(rep('ABC',1e6),collapse = ",")
  length(unlist(strsplit(cf$y,",")))
  for(i in 1:100){
    cat(i,"\n")
    h <- system.time(
      length(unlist(strsplit(cf$y,",")))
    )
    pl <- rbind(pl,
                data.frame(category=names(h)[1:3], tm=as.vector(unname(unlist(h)))[1:3],
                           it=i,elems='1M')
          )
    rm(h)
  }
  
  cf$y = paste0(rep('ABC',2e5),collapse = ",")
  length(unlist(strsplit(cf$y,",")))
  for(i in 1:100){
    cat(i,"\n")
    h <- system.time(
      length(unlist(strsplit(cf$y,",")))
    )
    pl <- rbind(pl,
                data.frame(category=names(h)[1:3], tm=as.vector(unname(unlist(h)))[1:3],
                           it=i,elems='200k')
    )
    rm(h)
  }
  
  cf$y = paste0(rep('ABC',5e6),collapse = ",")
  length(unlist(strsplit(cf$y,",")))
  for(i in 1:100){
    cat(i,"\n")
    h <- system.time(
      length(unlist(strsplit(cf$y,",")))
    )
    pl <- rbind(pl,
                data.frame(category=names(h)[1:3], tm=as.vector(unname(unlist(h)))[1:3],
                           it=i,elems='5M')
    )
    rm(h)
  }
  
  pl$elems <- paste0("n=",pl$elems)
  pl$tm <- pl$tm*1000
  pl$elems <- factor(pl$elems,levels=c("200k","1M", "5M"))
  
  qq<-ggplot(pl,aes(x=elems,y=log10(tm+1),fill=as.factor(elems)))+geom_violin()+
    geom_boxplot(width=0.1, fill="white",position = position_dodge(0.5))+
    facet_wrap(~category,nrow = 1,scales = 'free') +theme_bw()+
    theme(axis.text.x = element_text(size=14,colour = 'black',angle = 45,hjust = 1),
          axis.text.y = element_text(size=14,colour = 'black'),
          axis.title = element_text(size=15,colour = 'black'),
          legend.position = 'none',
          strip.text = element_text(size=15,colour = 'black'),
          strip.background = element_blank())+
    scale_y_continuous(breaks = c(0,1,2,3,4),labels = c(0,10,100,1000,10000))+
    xlab("")+ylab("Time \n(milliseconds)")
  pdf("C:/Dhawal/scRNASeq_data/Split_ops_Times.pdf",height = 4,width = 8)  
  qq
  dev.off()
  save(pl,file="C:/Dhawal/scRNASeq_data/Split_ops_Times.RData")
  
  
}

####------------------------------------------------------------------------------
####--------------- homo/ortho genes
####------------------------------------------------------------------------------
if(F){
  connGenes <- RSQLite::dbConnect(RSQLite::SQLite(),paste0("hs_genes.db"))
  RSQLite::dbListTables(connGenes)
  query <- paste0("SELECT * FROM gAliasProtein")
  gGenes <- RSQLite::dbGetQuery(connGenes, query)
  
  g <- read.delim("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt")
  g$DB.Class.Key
  g <- tidyr::separate_rows(g,Synonyms)
  g1 <- g[,-12] # - Synonyms
  g2 <- g[,-4] # - Symbol
  g2 <- cbind(g2[,c(1:3,11,4:10)])
  names(g2)[4] <- 'Symbol'
  names(g1)==names(g2)
  g1 <- unique(g1)
  g2 <- unique(g2)
  g <- rbind(g1,g2)
  rm(g1,g2)
  g <- unique(g)
  
  library(org.Hs.eg.db)
  columns(org.Hs.eg.db)
  keys(org.Hs.eg.db, keytype="ENSEMBL")[1:10]
  hs <- AnnotationDbi::select(org.Hs.eg.db, keys=keys(org.Hs.eg.db, keytype="ENSEMBL"),
                              keytype = "ENSEMBL",columns=c("SYMBOL","ENSEMBL"))
  g1 <- g[g$Common.Organism.Name=="human",]
  g1 <- g1[match(hs$SYMBOL,g1$Symbol),]
  g1$Symbol <- hs$ENSEMBL
  g1 <- g1[!is.na(g1$DB.Class.Key),]
  g <- rbind(g,g1)
  rm(g1,hs)

  library(org.Mm.eg.db)
  columns(org.Mm.eg.db)
  keys(org.Mm.eg.db, keytype="ENSEMBL")[1:10]
  mm <- AnnotationDbi::select(org.Mm.eg.db, keys=keys(org.Mm.eg.db, keytype="ENSEMBL"),
                              keytype = "ENSEMBL",columns=c("SYMBOL","ENSEMBL"))
  g1 <- g[g$Common.Organism.Name=="mouse, laboratory",]
  g1 <- g1[match(mm$SYMBOL,g1$Symbol),]
  g1$Symbol <- mm$ENSEMBL
  g1 <- g1[!is.na(g1$DB.Class.Key),]
  g <- rbind(g,g1)
  rm(g1,mm)
  
  library(org.Rn.eg.db)
  columns(org.Rn.eg.db)
  keys(org.Rn.eg.db, keytype="ENSEMBL")[1:10]
  rn <- AnnotationDbi::select(org.Rn.eg.db, keys=keys(org.Rn.eg.db, keytype="ENSEMBL"),
                              keytype = "ENSEMBL",columns=c("SYMBOL","ENSEMBL"))
  g1 <- g[g$Common.Organism.Name=="rat",]
  g1 <- g1[match(rn$SYMBOL,g1$Symbol),]
  g1$Symbol <- rn$ENSEMBL
  g1 <- g1[!is.na(g1$DB.Class.Key),]
  g <- rbind(g,g1)
  rm(g1,rn)
  
  g <- g[order(g$DB.Class.Key),]
  names(g)[4] <- 'geneSymbol'
  names(g)[7] <- "HGNC"
  names(g)[2] <- "Organism"
  g$species <- g$Organism
  g$species <- ifelse(g$Organism=="human","Homo sapien",g$species)
  g$species <- ifelse(g$Organism=="mouse, laboratory","Mus musculus",g$species)
  g$species <- ifelse(g$Organism=="rat","Rattus norvegicus",g$species)
  g$species <- ifelse(g$Organism=="zebrafish","Danio rerio",g$species)
  table(g$species)
  
  connGenes1 <- RSQLite::dbConnect(RSQLite::SQLite(),paste0("genes.db"))
  RSQLite::dbListTables(connGenes1)
  dbWriteTable(conn = connGenes1,name = "genehomologs",value = g)
  dbDisconnect(connGenes1)
  
  connGenes <- RSQLite::dbConnect(RSQLite::SQLite(),paste0("genes.db"))
  query <- paste0("SELECT * FROM genehomologs")
  gGenes <- RSQLite::dbGetQuery(connGenes, query)
  names(gGenes) <- c("DBClassKey", "Organism", "NCBITaxonID", "geneSymbol", 
                     "EntrezGeneID", "MGIID", "HGNC", "OMIMGeneID", "GeneticLocation", 
                     "GenomicCoordinates", "Name", "species") 
  gGenes$species <- gsub("Homo sapien","Homo sapiens",gGenes$species)
  table(gGenes$species)
  dbWriteTable(conn = connGenes,name = "genehomologs",value = gGenes,overwrite=T)
  dbDisconnect(connGenes)
  
  ##--- add gAliasProteins to genes.db for the time being
  connGenes <- RSQLite::dbConnect(RSQLite::SQLite(),paste0("hs_genes.db"))
  RSQLite::dbListTables(connGenes)
  query <- paste0("SELECT * FROM gAliasProtein")
  gGenes <- RSQLite::dbGetQuery(connGenes, query)
  connGenes1 <- RSQLite::dbConnect(RSQLite::SQLite(),paste0("genes.db"))
  RSQLite::dbListTables(connGenes1)
  dbWriteTable(conn = connGenes1,name = "gAliasProtein",value = gGenes)
  dbDisconnect(connGenes1)
  dbDisconnect(connGenes)
  
  
  g <- ifelse(g=="",NA,g)
  g$HGNC
  
  query <- paste0("SELECT * FROM gAliasProtein")
  gGenes <- RSQLite::dbGetQuery(connGenes, query)
  
  connGenes <- RSQLite::dbConnect(RSQLite::SQLite(),paste0("genes.db"))
  
  geneA = c("Itga5","Itgb1","a","dfgert")
  ORGANISM="Homo sapiens"
  gAlias <- local({
    gAlias <- list()
    for(f in geneA){
      gAlias[[f]] <- c(f)
      query <- paste0("SELECT * FROM genehomologs WHERE geneSymbol ='",f,"'")
      gGenes <- RSQLite::dbGetQuery(connGenes, query)$DBClassKey %>% unique()
      if(length(gGenes)==1){
        query <- paste0("SELECT * FROM genehomologs WHERE DBClassKey ='",gGenes,"'")
        gGenes <- RSQLite::dbGetQuery(connGenes, query)
        gGenes <- gGenes[gGenes$species==ORGANISM,]
        gAlias[[f]] <- unique(c(gAlias[[f]],gGenes$geneSymbol))
      }
    }
    gAlias
  })
  gAlias
  
  
}

####------------------------------------------------------------------------------
####--------------- using efficient data storage solution
####------------------------------------------------------------------------------
if(F){
  rm(list=ls())
  gc()
  setwd("/home/rstudio/scripts/SciViewerDev/SciViewer_V1/")
  DATADIR="/home/rstudio/data/SCS/"
  source("plot_functions.R")
  connList <- list(A=RSQLite::dbConnect(RSQLite::SQLite(),paste0(DATADIR,"scRNA_RNA_HS_healthyCholangiocytes_10x.db")))#,
                  # B=RSQLite::dbConnect(RSQLite::SQLite(),paste0(DATADIR,"scRNA_LungDB02.db")),
                  # C=RSQLite::dbConnect(RSQLite::SQLite(),paste0(DATADIR,"scRNA_Heart01.db")))
  connGenes <- RSQLite::dbConnect(RSQLite::SQLite(),paste0("genes.db"))
  
  
  ##------ update studies
  sc_studies <- c()
  for(i in 1:length(connList)){
    tablist <- RSQLite::dbListTables(connList[[i]])
    tablist <- tablist[grep("_study",tablist)]
    for(j in tablist){
      cat(i,"\t",j,"\n")
      query <- paste0("SELECT * FROM ",j)
      sc_studies <- rbind(sc_studies,
                          cbind(RSQLite::dbGetQuery(connList[[i]], query), 
                                ObjID=names(connList)[i]))
    }
    rm(i,j)
  }
  
  #--- avail genes
  sc_study <- sc_studies$Database[1]
  connID <- sc_studies[sc_studies$Database==sc_study,]$ObjID %>% as.character
  query <- paste0("SELECT geneSymbol FROM ",sc_study,"_data")
  cat(" getting gene lists for: ", sc_study, " objID: ", connID,"\n")
  availGenes <- queryDB(HANDLER=connList[[connID]], 
                        QUERY=query,REPO_NAME=NULL,
                        USE_REMOTE_DB=FALSE)
  availGenes <- tryCatch({
    query <- paste0("SELECT * FROM genehomologs WHERE species = '",sc_studies[sc_studies$Database==sc_study,]$ORGANISM,"'")
    g1 <- RSQLite::dbGetQuery(connGenes, query)
    g1 <- g1[,c("geneSymbol","DBClassKey")] %>% unique()
    availGenes <- merge(availGenes,g1,by="geneSymbol",all.x=T)
    g1 <- g1[-grep("^ENS",g1$geneSymbol),]
    names(g1)[1] <- "Name"
    availGenes <- merge(availGenes,g1,by="DBClassKey",all.x=T)
    availGenes$Name <- ifelse(is.na(availGenes$Name),as.character(availGenes$geneSymbol),as.character(availGenes$Name))
    availGenes
  },error=function(e){
    availGenes <- cbind(DBClassKey=NA,availGenes,Name=availGenes$geneSymbol)
    availGenes
  })
  #availGenes <- as.character(unique(availGenes$geneSymbol))
  geneA = availGenes$geneSymbol[1]
  
  sum(is.na(availGenes$Name))

  ##------ update genes
  #gAlias <- local({
  #  query <- paste0("SELECT * FROM gAliasProtein")
  #  gGenes <- RSQLite::dbGetQuery(connGenes, query)
  #  gAlias <- list()
  #  for(f in geneA){
  #    hgnc <- gGenes[gGenes$geneSymbol==f,]$HGNC %>% as.character
  #    gg <- as.character(gGenes[gGenes$HGNC%in%hgnc,]$geneSymbol)
  #    gAlias[[f]] <- unique(c(f,gg))
  #  }
  #  gAlias
  #})
  gAlias <- local({
    gAlias <- list()
    gg <- availGenes
    for(f in geneA){
      #gAlias[[f]] <- unique(c(f,gg[gg$Name==f,]$geneSymbol[1]))
      gAlias[[f]] <- f
    }
    gAlias
  })
  cat("Alias gene names: ", paste0(unlist(gAlias),collapse = ";"),"\n")
  
  
  ##--- louvain
  query <- paste0("SELECT * FROM ",sc_study,"_metaFeatures")
  louvain <- RSQLite::dbGetQuery(connList[[connID]], query)
  louvain$value <- NA
  
  ##-- avail features
  feature <- features <- unlist(strsplit(sc_studies$CATVAR[1],","))[1]
  add_features <- unlist(strsplit(sc_studies$CATVAR[1],","))
  add_features <- add_features[add_features!='SVcell_type']
  genename <-  availGenes$geneSymbol[1]
  
  
  
  #--- summary count plots
  dhc_columnPlot(pl=louvain,
                 features = add_features[2],
                 ycol = 'value',plotType = 'count',
                 log.transform = F,main = paste0("Cell type proportions<br>(",add_features[2],")"),
                 xlab = "",ylab = "% of total",
                 sourceref = 'SciViewer')
  
  
  
  ## cross study comparison
  mygene=geneA
  for(i in 1:nrow(sc_studies)){
    pl <- NULL
    ql <- NULL
    gl <- NULL
    cntr =1 
    cf <- sc_studies[i,]
    cat(" extracting data for the gene",mygene," in ",cf$Database,"\n")
    while(cntr<=length(gAlias[[mygene]]) && !isTRUE(nrow(pl)>0)){
      cat("  gene/alias name: ", gAlias[[mygene]][cntr],"\n")
      ##-- 01) Expression
      pl <- local({
        query <- paste0("SELECT * FROM ",as.character(cf$Database),"_FeatureSummary WHERE feature = '",mygene,"'")
        tryCatch({
          queryDB(HANDLER=connList[[as.character(cf$ObjID)]], 
                  QUERY=query)
        },error=function(e){ NULL
        })
      })
      ##-- 02) Marker
      gl <- local({
        query <- paste0("SELECT * FROM ",as.character(cf$Database),"_Marker WHERE geneSymbol = '",mygene,"'")
        tryCatch({
          queryDB(HANDLER=connList[[as.character(cf$ObjID)]], 
                  QUERY=query)
        },error=function(e){ NULL
        })
      })
      ##-- 03) BioMarker
      ql <- local({
        query <- paste0("SELECT * FROM ",as.character(cf$Database),"_DEG WHERE geneSymbol = '",mygene,"'")
        tryCatch({
          queryDB(HANDLER=connList[[as.character(cf$ObjID)]], 
                  QUERY=query)
        },error=function(e){ NULL
        })
      })
      ## counter
      cntr = cntr +1 
    }
    if(isTRUE(nrow(pl)>0)){
      pl$Database <- as.character(cf$Database)
      #compvariables$geneExp <- rbind(compvariables$geneExp, pl)
      cat("Expression Comparison: A study added\n")
    }
    if(isTRUE(nrow(gl)>0)){
      gl$Database <- as.character(cf$Database)
      #compvariables$geneMarker <- rbind(compvariables$geneMarker, gl)
      cat("Marker Expression Comparison: A study added\n")
    }
    if(isTRUE(nrow(ql)>0)){
      ql$Database <- as.character(cf$Database)
      #compvariables$geneBioMarker <- rbind(compvariables$geneBioMarker, ql)
      cat("Biological Marker Expression Comparison: A study added\n")
    }
  }
  

  ##--- gene expression dataframes
  pldf<-list()
  genenames = availGenes$geneSymbol[1:5]
  for(genename in genenames){
    pldf[[genename]] <- get_plot_df_sc04(connSc = connList[[connID]],
                                                  sc_study,
                                                  genename=genename,
                                                  louvain,REPO_NAME=NULL,USE_REMOTE_DB=FALSE)  
    #pldf[[genename]]$value <- ifelse(pldf[[genename]]$value<0,0,round(pldf[[genename]]$value))
  }
  #qq <- plot_multigene_grouped_heatmap(pldf = pldf,genenames = genename,feature)
  #qq[2]
  
  plot_umap_raster(pl=pldf[[1]],
                        feature = "cell_type",
                        genename = "A",
                        legendPlot = F,
                        selectedFeature = NULL,
                        source.val=NULL,js=NULL)
  
  ## umap
  js <- "
    function(el, x, inputName){
      var id = el.getAttribute('id');
      var gd = document.getElementById(id);
      var xScrollPos = el.scrollLeft || document.documentElement.scrollLeft;
      var yScrollPos = el.scrollTop || document.documentElement.scrollTop;
      var d3 = Plotly.d3;
      console.log(Plotly);
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
  p <- plot_umap_raster(pl=pldf[[1]],
                   feature = "cell_type",
                   genename = "A",
                   legendPlot = F,
                   selectedFeature = NULL,
                   source.val=NULL,js=NULL)
  
  p <- plotly::config(p,mathjax = 'cdn') %>% partial_bundle(local = F)
  
  #withMathJax()
  #p[["dependencies"]][[5]][["script"]] <- "https://cdn.plot.ly/plotly-2.16.1.min.js"
  p %>% htmlwidgets::onRender(jsCode = js, 
                              data = environment(p[["x"]][["visdat"]][[1]])[["source"]])
  
  p$dependencies[[5]]$script = "https://cdn.plot.ly/plotly-latest.min.js"
  
  system.file(package = "plotly", "htmlwidgets/lib/plotlyjs")
  p$dependencies[[5]]$src$file = "https://cdn.plot.ly/plotly-latest.min.js"
  htmltools::htmlDependency()
  
  p %>% 
    htmlwidgets::onRender("
    function(el) { 
      el.on('plotly_hover', function(d) { 
        console.log('Hover: ', d); 
      });
      el.on('plotly_click', function(d) { 
        console.log('Click: ', d);
      });
      el.on('plotly_selected', function(d) { 
        console.log('Select: ', d); 
      });
    }
  ")
  
  q <- plot_ly(mtcars, x = ~wt, y = ~mpg) 
  q %>% layout(showLegend=T)
  
  q %>% 
    htmlwidgets::onRender("
    function(el) { 
      el.on('plotly_hover', function(d) { 
        console.log('Hover: ', d); 
      });
      el.on('plotly_click', function(d) { 
        console.log('Click: ', d);
      });
      el.on('plotly_selected', function(d) { 
        console.log('Select: ', d); 
      });
    }
  ")
  

  ##------ scatter
  plot_multi_umap_raster(list(pldf[[1]]),
                         genenames = "a",
                         legendPlot = F,xrange = NULL,yrange = NULL)
  
  ##--- correlation scatter
  cf <- data.frame(pldf[[1]][42],pldf[[2]][42])
  names(cf) <- names(pldf)
  cf$variable <- apply(cf,1,function(x) ifelse(sum(is.na(x))==0,1,0) )
  cf <- as.matrix(cf)
  cf <- ifelse(is.na(cf),0,cf)
  cf <- as.data.frame(cf)
  cf <- cbind(celltype=pldf[[1]]$cell_type,cf)
  cf$value <- cf$variable
  dhc_columnPlot(pl = cf,features = c("celltype","variable"),plotType = "count",ycol = "value")
  
  pl <- louvain
  pl$value <- cf$value
  
  plot_multi_umap_raster(list(pl),
                         genenames = "co-expressing cells",
                         legendPlot = F,xrange = NULL,yrange = NULL)
  

  smoothScatter(cf$x,cf$y)
  #plotRasterly(cf,mapping = aes(x=x,y=y),drop_data = T,as_image = T)

  
  
  
  
  # summarize by donor
  pl <- pldf[[1]]
  pl <- pl[!is.na(pl$value),]
  dhc_columnPlot(pl=pl,
                 features = 'cell_type',
                 plotType = 'expn',min.cell.number = 10,
                 ycol = 'value',
                 log.transform = F,main = paste0(" expression"),
                 xlab = "",ylab = "Normalized expression",
                 sourceref= "SOURCEREF")
  
  
  feature = c('cell_type',"donor_uuid")
  ycol='value'
  gl <- pldf[[1]] %>% dplyr::group_by_at(feature) %>% dplyr::select(dplyr::all_of(ycol)) %>% 
    dplyr::summarise_all(c("avg_expn"=mean))

  dhc_columnPlot(pl=gl,
                 features = 'cell_type',
                 plotType = 'expn',min.cell.number = 0,
                 ycol = 'avg_expn',
                 log.transform = F,main = " expression",
                 xlab = "",ylab = "Normalized expression",
                 sourceref="SOURCEREF")
  
  
  
  
  ##--- markers
  query <- paste0("SELECT * FROM ",sc_study,"_Marker  WHERE geneSymbol = '","ENSG00000171612","'")
  pl <- tryCatch({
    queryDB(HANDLER=connList[[connID]], 
            QUERY=query)
  },error=function(e){
    return(NULL)
  })
  
  RSQLite::dbGetQuery(connList[[connID]], query)
  
  
  ##------ barplot
  conqpl <- c()
  for(i in 1:length(pldf)){
    conqpl <- rbind(conqpl,cbind(gName=names(pldf)[[i]],
                                 pldf[[i]][,c('cell_type',"value")]
                                 )
                    )
  }
  dhc_columnPlot(pldf[[1]], features=c('SVcell_type'),ycol='value',
                 log.transform = F,
                 main="",xlab="",
                 ylab="average scaled expression <br/>(log2)",
                 min.cell.number=0,
                 plotType='expn')
  dhc_boxPlot(pldf[[1]], features=c('cell_type',add_features[2]),ycol='value',
    log.transform = F,
    main="",xlab="",
    ylab="average scaled expression <br/>(log2)"
  )
  
  ##--- dot plot
  
  
  ##--- corerelation
  z <- data.frame(cell_type=pldf[[1]]$cell_type)
  for(i in 1:length(pldf)){
    cat(" corr matrix: i=",i," gene=",names(pldf)[i],"\n")
    z <- cbind(z,pldf[[i]][,'value'])
  }
  names(z)[2:ncol(z)] <- names(pldf)
  z <- z[complete.cases(z),]
  z <- z[z$cell_type=='Lymphatic',]
  z <- z[,2:ncol(z)]
  flt <- apply(z, 1, function(x){
    ifelse(sum(x==0)==length(x),F,T)
  })
  z <- z[flt,]
  cat("ssc: corr matrix dataframe dim: ", nrow(cor(z)),"\n")
  hchart_cor(cor(z))
  
  ##---- scaling
  colScale <- c("#fcfbfd", "#9e9ac8", "#3f007d")
  minValCol <- "#f0f0f0"
  x='V1'
  y='V2'
  ycol='value'
  val.range <- c()
  nuniques <- c()
  uqvals <- c()
  for(i in 1:length(pldf)){
    pldf[[i]][,ycol] <- as.numeric(pldf[[i]][,ycol])
    val.range <- c(val.range,range(pldf[[i]][,ycol],na.rm=T))
    uqvals <- c(uqvals, unique(pldf[[i]][,ycol] ))
  }
  uqvals <- sort(unique(uqvals),decreasing = F)
  val.range <- range(val.range,na.rm = T)
  nuniques <- length(uqvals)
  pal <- grDevices::colorRampPalette(colScale)(nuniques)
  pal[1] <- minValCol
  coldf <- data.frame(val=uqvals,col=pal,stringsAsFactors = F)
  
  i=4
  gl <- pldf[[i]][,c(x,y,ycol)]
  names(gl) <- c("x","y","feature")
  xx <- unique(gl$feature)
  xx <- coldf[coldf$val%in%xx,]
  xx <- xx[order(xx$val,decreasing = F),]
  
  #xx$val <- paste0("x",xx$val)
  #gl$feature <- paste0("x",gl$feature)
  #gl$feature <- factor(gl$feature,levels=xx$val)
  #ggRasterly(gl,aes(x=x,y=y,color=feature))
  #ggplot(gl[1:100000,],aes(x,y,col=feature))+geom_point()
  
  gl[1:100,]%>%
    plotRasterly(aes(x=x,y=y,on=feature),
                 color=colScale,reduction_func='sum',
                 show_raster = T,drop_data = T,background = 'white',
                 legend=T,
                 as_image = T,variable_check = F,
                 plot_width = 300,plot_height = 300) %>% suppressMessages()
  
  
  gl[1:100,]%>%
    rasterly(aes(x=x,y=y,on=feature),color=(xx$col),
           show_raster = T,drop_data = T,variable_check = T)%>% 
    rasterly_points() %>% 
  rasterly_build()
  
  color=c("0"="#f0f0f0", "1"="#F0EEF6", "2"="#E4E2EF", "3"="#D8D6E9", "4"="#CDCAE2", "5"="#C1BEDB",
          "6"="#B5B2D5", "7"="#A9A6CE", "8"="#9E9AC8", "9"="#9286BE", "10"="#8673B5", "11"="#7A60AB", 
          "12"="#6E4DA2", "13"="#623999", "14"="#56268F", "15"="#4A1386", "16"="#3F007D")  
  
  
  range(pl$value)
  hist(pl$value,breaks = 100)
  
  #pl$Donor <- as.character(pl$Donor)
  #pl$Donor <- paste0("donor",pl$Donor)
  dhc_columnPlot(pl,features = c('Donor','cell_group'),ycol = 'value',
                 plotType = 'expn',min.cell.number = 10,
                 log.transform = F,main = paste0(genename," expression"),
                 xlab = "",ylab = "% of total")
  
  
  
  
  d <- read.delim("C:/Dhawal/scRNASeq_data/Kaminski_IPF2020/GSE136831_AllCells.cellBarcodes.txt.gz")
  
}

####------------------------------------------------------------------------------
####--------------- trash
####------------------------------------------------------------------------------
if(F){
  rm(list=ls())
  setwd("C:/Dhawal/scripts/pDDSingleCellApp")
  DATADIR="C:/Dhawal/SHINYAPP_DATA/"
  source("plot_functions.R")
  connList <- list(#A=RSQLite::dbConnect(RSQLite::SQLite(),paste0(DATADIR,"scRNA_ChuaNat2020.db")))#,
   B=RSQLite::dbConnect(RSQLite::SQLite(),paste0(DATADIR,"scRNA_LungDB03.db")),
   C=RSQLite::dbConnect(RSQLite::SQLite(),paste0(DATADIR,"scRNA_LungDB02.db")))
  connGenes <- RSQLite::dbConnect(RSQLite::SQLite(),paste0(DATADIR,"hs_genes.db"))
  
  
  ##------ update genes
  query <- paste0("SELECT * FROM gAliasProtein")
  gGenes <- RSQLite::dbGetQuery(connGenes, query)
  availGenes <- as.character(unique(gGenes$geneSymbol))
  
  
  
  ##------ update studies
  sc_studies <- c()
  for(i in 1:length(connList)){
    tablist <- RSQLite::dbListTables(connList[[i]])
    tablist <- tablist[grep("_study",tablist)]
    for(j in tablist){
      cat(i,"\t",j,"\n")
      query <- paste0("SELECT * FROM ",j)
      sc_studies <- rbind(sc_studies,
                          cbind(RSQLite::dbGetQuery(connList[[i]], query), 
                                ObjID=names(connList)[i]))
    }
    rm(i,j)
  }
  
  study <- sc_studies$Database[1]
  connID <- sc_studies[sc_studies$Database==study,]$ObjID %>% as.character
  feature <- features <- unlist(strsplit(sc_studies$CATVAR[1],","))[1]
  add_features <- unlist(strsplit(sc_studies$CATVAR[1],","))
  add_features <- add_features[add_features!='cell_type']
  genename <-  "YY1"
  
  query <- paste0("SELECT * FROM ",study,"_DGE")
  pl <-  RSQLite::dbGetQuery(connList[[connID]], query)
  head(pl)
  #c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Test", 
  #  "Tag", "geneSymbol")
  
  
  connSc <- RSQLite::dbConnect(RSQLite::SQLite(), paste0("scRNA_Heart01.db"))
  dbListTables(connSc)
  
  study="AA"
  d <- read.delim("C:/Dhawal/scRNASeq_data/PCL/diffexpress_human_ICM_vs_NF_differential_expression_v2.0.txt",header=T)
  d$Test <- paste0(d$grp1,"-",d$grp2)
  d$B <- 0
  d <- d[,c("logFC_CB", 
            "AveExpr_CB", "t_CB", "P.Value_CB", "adj.P.Val_CB","B","Test",
            "celltype","gene", "ensembl_id","grp1", "grp2", "nsamp_grp1", 
            "nsamp_grp2", "ncell_grp1", "ncell_grp2", "pctgt0_grp1", "pctgt0_grp2", 
            "avgexpr_grp1", "avgexpr_grp2", "gene_bkg_prob", "logFC_CR", 
            "AveExpr_CR", "t_CR", "P.Value_CR", "adj.P.Val_CR")]
  names(d)[1:9] <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Test", 
                       "Tag", "geneSymbol")
  RSQLite::dbWriteTable(connSc,"_DGE",gGenes)
  dbDisconnect(connSc)
  
  connSc <- RSQLite::dbConnect(RSQLite::SQLite(), paste0("scRNA_Heart01.db"))
  dbListTables(connSc)
  study="AA"
  d <- read.delim("C:/Dhawal/scRNASeq_data/PCL/human_dcm_hcm_freeze3_disease_vs_healthy_de_hcm_vs_nf_de.txt",header=T)
  e <- read.delim("C:/Dhawal/scRNASeq_data/PCL/human_dcm_hcm_freeze3_disease_vs_healthy_de_dcm_vs_nf_de.txt",header=T)
  names(d)==names(e)
  d <- rbind(d,e)
  rm(e)
  d$Test <- paste0(d$grp1,"-",d$grp2)
  d$B <- 0
  d <- d[,c("logFC_CB", 
            "AveExpr_CB", "t_CB", "P.Value_CB", "adj.P.Val_CB","B","Test",
            "celltype","gene", "ensembl_id","grp1", "grp2", "nsamp_grp1", 
            "nsamp_grp2", "ncell_grp1", "ncell_grp2", "pctgt0_grp1", "pctgt0_grp2", 
            "avgexpr_grp1", "avgexpr_grp2", "flag_bkg", "flag_lowexpress",
            "logFC_CR", 
            "AveExpr_CR", "t_CR", "P.Value_CR", "adj.P.Val_CR")]
  names(d)[1:9] <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Test", 
                     "Tag", "geneSymbol")
  
  query <- paste0("SELECT * FROM HCM_metaFeatures")
  meta <-  RSQLite::dbGetQuery(connSc, query)
  dput(names(table(d$Tag)))
  d$Tag <- ifelse()
  
  c("cardiomyocyte","fibroblast","neuronal", "lymphatic_endothelial",
    "endocardial", "endothelial1", 
    "endothelial2","adipocyte",  "endothelial3",   
    "lymphocyte1", "lymphocyte2", "macrophage",  "pericyte", 
    "pseudo_bulk", "vsmc")
  c("0: Cardiomyocyte", "1: Fibroblast", "10: Neuronal", "11: Lymphatic Endothelial",
    "12: Cardiomyocyte II", "13: Fibroblast II", "14: Unknown/Garbage",
    "2: Endothelial I", "3: Mural Cells", "4: Macrophage", "5: Lymphocyte I",
    "6: Endothelial II", "7: Adipocyte", "8: Activated Fibroblast",
    "9: Lymphocyte II")
  
  
  RSQLite::dbWriteTable(connSc,"_DGE",gGenes)
  dbDisconnect(connSc)
  
}



