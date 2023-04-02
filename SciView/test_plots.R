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
  ## JAX homolog database
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
  
  ##-- HomologGene database
  m <- read.delim("https://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data",header=F)
  #--- manually prepared db
  taxdf <- data.frame(taxid=c(318829,33169,28985,10116,10090,9913,9615,9606,9598,9544,9031,8364,7955,7227,7165,6239,5141,4932,4896,4530,3702),
                      taxname=c("Pyricularia oryzae","Eremothecium gossypii","Kluyveromyces lactis","Rattus norvegicus","Mus musculus","Bos taurus","Canis lupus familiaris","Homo sapiens","Pan troglodytes","Macaca mulatta","Gallus gallus","Xenopus tropicalis","Danio rerio","Drosophila melanogaster","Anopheles gambiae","Caenorhabditis elegans","Neurospora crassa","Saccharomyces cerevisiae","Schizosaccharomyces pombe","Oryza sativa","Arabidopsis thaliana"))
  taxdf <- taxdf[match(m$V2,taxdf$taxid),]
  m$Common.Organism.Name <- taxdf$taxname
  rm(taxdf)
  names(m) <- c("V1","NCBI.Taxon.ID","EntrezGene.ID","Symbol","V5","V6","Common.Organism.Name")
  gx <- g[,c("DB.Class.Key","EntrezGene.ID")] %>% unique()
  gx <- gx[match(m$EntrezGene.ID,gx$EntrezGene.ID),]
  m$DB.Class.Key <- gx$DB.Class.Key
  mx <- m[,c("V1","DB.Class.Key")] %>% unique()
  mx <- mx[complete.cases(mx),]
  mx <- mx[match(m$V1,mx$V1),]
  m$DB.Class.Key <- mx$DB.Class.Key
  mx <- m[,c("V1","DB.Class.Key")] %>% unique()
  my <- mx[is.na(mx$DB.Class.Key),]
  mx <- mx[complete.cases(mx),]
  range(mx$DB.Class.Key)
  my$DB.Class.Key <- sample(45000000:50000000,length(my$V1),replace = F)
  length(unique(my$DB.Class.Key))
  mx <- rbind(mx,my)
  mx <- mx[match(m$V1,mx$V1),]
  m$DB.Class.Key <- mx$DB.Class.Key
  m <- m[,c("DB.Class.Key","Common.Organism.Name", "NCBI.Taxon.ID",
            "Symbol","EntrezGene.ID")]
  m$Mouse.MGI.ID <- NA
  m$HGNC.ID <- NA
  m$OMIM.Gene.ID <- NA
  m$Genetic.Location <- NA
  m$Genomic.Coordinates..mouse....human... <- NA
  m$Name <- NA
  names(m)==names(g)
  m$tmp <- paste0(m$DB.Class.Key,m$NCBI.Taxon.ID)
  m <- m[!m$tmp%in%paste0(g$DB.Class.Key,g$NCBI.Taxon.ID),]
  m$tmp = NULL
  g <- rbind(g,m)
  g <- unique(g)
  rm(m,mx,my,gx,)
    
  library(org.Hs.eg.db)
  columns(org.Hs.eg.db)
  keys(org.Hs.eg.db, keytype="ENSEMBL")[1:10]
  hs <- AnnotationDbi::select(org.Hs.eg.db, keys=keys(org.Hs.eg.db, keytype="ENSEMBL"),
                              keytype = "ENSEMBL",columns=c("SYMBOL","ENSEMBL"))
  g1 <- g[g$NCBI.Taxon.ID==9606,]
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
  g1 <- g[g$NCBI.Taxon.ID==10090,]
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
  g1 <- g[g$NCBI.Taxon.ID==10116,]
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
  g$species <- ifelse(g$Organism=="Homo sapiens","Homo sapien",g$species)
  table(g$species)
  
  
  ##--- write database file
  require(RSQLite)
  connGenes <- RSQLite::dbConnect(RSQLite::SQLite(),paste0("genes.db"))
  RSQLite::dbListTables(connGenes)
  dbWriteTable(conn = connGenes,name = "genehomologs",value = g,overwrite=T)
  dbDisconnect(connGenes)
  
  
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
  setwd("/home/rstudio/scripts/SciViewerDev/SciView/")
  DATADIR="/home/rstudio/data/SCS/final/"
  source("plot_functions.R")
  source("vars.R")
  connList <- list(A=RSQLite::dbConnect(RSQLite::SQLite(),paste0(DATADIR,"scRNA_RNA_HS_DCMAdipocytes_10x.db")),
                   B=RSQLite::dbConnect(RSQLite::SQLite(),paste0(DATADIR,"scRNA_RNA_HS_IlleumChrohnsDisease_10X.db")),
                   C=RSQLite::dbConnect(RSQLite::SQLite(),paste0(DATADIR,"scRNA_RNA_Hs_TabulaSapiensEpithelial_10x.db")),
                   D=RSQLite::dbConnect(RSQLite::SQLite(),paste0(DATADIR,"scRNA_RNA_HS_CrohnDiseaseCCF_10x.db"))
                   )
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
    availGenes <- availGenes[,c(2,3,1)]
    availGenes
  },error=function(e){
    availGenes <- cbind(DBClassKey=NA,availGenes,Name=availGenes$geneSymbol)
    availGenes
  })
  #availGenes <- as.character(unique(availGenes$geneSymbol))
  #geneA = c('ITGB8','ITGAV','CFBP')##availGenes$geneSymbol[6]
  geneA = availGenes$geneSymbol[6]
  
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
      gAlias[[f]] <- unique(c(f,gg[gg$Name==f,]$geneSymbol[1]))
      #gAlias[[f]] <- f
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
  
  
  ##--- DEG marker
  pl <- tryCatch({
    query <- paste0("SELECT * FROM ",sc_studies$Database[1],"_DEG")
    queryDB(HANDLER=connList[[connID]], 
            QUERY=query,REPO_NAME=REPO_NAME,
            USE_REMOTE_DB=USE_REMOTE_DB)
  },error=function(e){
    NULL
  })
  
  
  ##--- gene as marker
  if(F){
    cf=sc_studies[1,]
    gf <- local({
      gf <- c()
      for(g in geneA){
        pl <- c()
        cntr=1
        cat("querying gene for celltype marker: ",g,"\n")
        while(cntr<=length(gAlias[[g]]) && !isTRUE(nrow(pl)>0)){
          pl <- tryCatch({
            query <- paste0("SELECT * FROM ",cf$Database,"_Marker WHERE geneSymbol = '",gAlias[[g]][cntr],"'")
            queryDB(HANDLER=connList[[connID]], 
                    QUERY=query,REPO_NAME=REPO_NAME,
                    USE_REMOTE_DB=USE_REMOTE_DB)
          },error=function(e){
            NULL
          })
          cntr = cntr + 1
        }
        if(isTRUE(nrow(pl)>0)){
          gf <- rbind(gf,pl)
        }
        rm(pl)
      }
      gf
    })
  }
  
  
  
  ##--- cell counts/proportions
  if(F){
    z <- louvain %>% 
      dplyr::group_by_at(add_features[c(2,7,8)]) %>% 
      dplyr::select(dplyr::all_of('value')) %>% 
      dplyr::summarise_all(c(length))
    names(z) <- c("group","subgroup","donor" ,"count")
    z <- z %>%
      group_by(donor) %>%
      mutate(subgroupcount = sum(count))
    z$fraction <- round(z$count/z$subgroupcount,3) 
    
    dhc_boxPlot(pl=z,
                features = c('group','subgroup'),
                ycol = 'fraction',
                log.transform = F,
                main = paste0(" fraction"),
                xlab = "",ylab = 'ylab',
                sourceref='SOURCEREF')
    
  }
  
  
  
  ## compare expression across studies
  USE_REMOTE_DB = FALSE
  REPO_NAME="Minuteman"
  pl <-c()
  for(i in 1:3){
    cf<-sc_studies[i,]
    query <- paste0("SELECT * FROM ",as.character(cf$Database),"_FeatureSummaryKeys")
    x <- queryDB(HANDLER=connList[[as.character(cf$ObjID)]], 
                 QUERY=query,REPO_NAME=REPO_NAME,
                 USE_REMOTE_DB=USE_REMOTE_DB)
    query <- paste0("SELECT * FROM ",as.character(cf$Database),"_FeatureSummary WHERE feature = '",geneA,"'")
    y <- queryDB(HANDLER=connList[[as.character(cf$ObjID)]], 
                 QUERY=query,REPO_NAME=REPO_NAME,
                 USE_REMOTE_DB=USE_REMOTE_DB)
    x$norm_avg_priorLT <- unname(unlist(y[,2:ncol(y)]))
    x$feature <- unname(unlist(unique(y[1,1])))
    message(paste0(names(x),collapse = ","))
    x$Database <- as.character(cf$Database)
    pl <- rbind(pl,x)
    rm(y,x,i,query,cf)
  }
  
  
  myfun01 <- function(x){ log2(median(2^x,na.rm=T)) }
  qf <- pl %>%
    dplyr::group_by(cell_type,feature,Database) %>%
    dplyr::summarize(MedExpr=myfun01(norm_avg_priorLT),expression_values = list(norm_avg_priorLT)) %>%
    mutate(boxplot = NA)
  qf$expression_values <- lapply(qf$expression_values, function(i){
    x <- ifelse(length(i)==1,paste0(rep(i,2),collapse = " "),paste0(i,collapse = " "))
    as.numeric(unlist(strsplit(x," ")))
  })
  reactable(qf, 
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
  

  ##------ scatter
  plot_multi_umap_raster(list(pldf[[1]]),
                         genenames = "a",
                         legendPlot = F,xrange = NULL,yrange = NULL)
  
  ##--- ridgeplot
  z <- data.frame(pldf[[1]][40],pldf[[2]][40])
  names(z) <- names(pldf)
  z <- as.matrix(z)
  z[is.na(z)] <- 0
  z <- as.data.frame(z)
  z <- cbind(celltype=pldf[[1]]$cell_type,z)
  mx <- reshape::melt(z,measure.vars=names(z[,-1]))
  ggplot(mx, aes(x = value, y = variable,fill=stat(x))) +
    geom_density_ridges_gradient(scale = 3,alpha=0.5) +
    scale_fill_viridis_c(name = "Avg. Expn",space = 'Lab',alpha = .3,option = 'C') +
    labs(title = 'Gene expression')+xlab("Expression value")+ylab("")+
    theme_bw()+
    theme(axis.text = element_text(color = 'black'),
          legend.position = 'bottom')
  
  
  ##--- 2ndary gene expression across 
  BASEVAL=0
  GENE0 = "ENSG00000181924"
  z <- data.frame(pldf[[1]][40],pldf[[2]][40])
  names(z) <- names(pldf)
  z <- as.matrix(z)
  z[is.na(z)] <- 0
  z <- as.data.frame(z)
  z <- cbind(celltype=pldf[[1]]$cell_type,z)
  z$basegene <- ifelse(z[,GENE0]>BASEVAL,"Expressing","Non-expressed")
  z[,GENE0] <- NULL
  z <- reshape::melt(z,measure.vars=names(z)[!names(z)%in%c('celltype','basegene')])
  z %>%
    ggplot(aes(x=celltype,y=value,fill=basegene))+
    geom_violin(scale = 'width',draw_quantiles = c(0.5))+
    facet_wrap(~variable,ncol = 1,scales = "free_y",strip.position = "right")+
    theme_bw()+theme(legend.position = "right",
                     axis.text = element_text(size=14),
                     axis.text.x = element_text(angle=-45,hjust=0,vjust=1),
                     axis.title = element_text(size=16),
                     strip.text = element_text(angle=0,hjust=1,size = 14),
                     strip.background = element_blank(),
                     plot.title = element_text(size=20,hjust=0.5),
                     plot.subtitle = element_text(size=12,hjust=0.5))+
    xlab("")+ylab("norm.expresion")+
    guides(fill=guide_legend(title=GENE0))
  
  
  
  ##-- column plot
  z$value <- z$variable
  dhc_columnPlot(pl = cf,features = c("celltype","variable"),plotType = "count",ycol = "value")
  
  
  
  ##--- correlation scatter
  cf <- data.frame(pldf[[1]][40],pldf[[2]][40])
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

if(F){
  library(dplyr)
  library(sparkline)
  library(reactable)
  
  data <- chickwts %>%
    group_by(feed) %>%
    summarise(weight = list(weight)) %>%
    mutate(boxplot = NA, sparkline = NA)
  
  data$weight[[6]] <- c(123)
  data$weight[[6]] <- c(123,NA)
  
  
  data$weight <- lapply(data$weight, function(i){
    x <- ifelse(length(i)==1,paste0(rep(i,2),collapse = " "),paste0(i,collapse = " "))
    as.numeric(unlist(strsplit(x," ")))
  })
  
  
  
  reactable(data, columns = list(
    feed = colDef(filterable = T,searchable = T,sortable = T,resizable = T),
    weight = colDef(cell = function(values) {
      sparkline(values, type = "bar", chartRangeMin = 0, chartRangeMax = max(chickwts$weight))
    }),
    boxplot = colDef(cell = function(value, index) {
      sparkline(data$weight[[index]], type = "box")
    },sortable = T,resizable = T),
    sparkline = colDef(cell = function(value, index) {
      sparkline(data$weight[[index]])
    },show = F)
  ))
}



if(F){
  
  so <- readInpscRNAseq(filepath="/home/rstudio/data/SCS/Hs_TabulaSapiensEpithelial_10x.h5ad",intractv = F)
  range(so[["RNA"]]@data)
  
  m <- AverageExpression(object = so,features = 'ENSG00000105855',return.seurat = FALSE,
                    group.by = c('cell_type','donor_id'),slot = "data",verbose = TRUE)
  m <- m[[1]] %>% as.data.frame()
}
