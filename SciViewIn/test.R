rm(list=ls())
setwd("~/scripts/SciViewerDev/SciViewIn/")
gc()

##--- code to reformat .gmt file; 
if(F){
  
  connGenes <- RSQLite::dbConnect(RSQLite::SQLite(),'C:/Dhawal/SHINYAPP_DATA/hs_genes.db')
  RSQLite::dbListTables(connGenes)
  query <- paste0("SELECT * FROM gAlias")
  genes <- RSQLite::dbGetQuery(connGenes, query)
  RSQLite::dbDisconnect(connGenes)
  rm(connGenes,query)
  
  out <- c()
  for(myfile in c("c2.cp.kegg.v7.5.1.symbols.gmt","c2.cp.wikipathways.v7.5.1.symbols.gmt")){
    d <- GSA::GSA.read.gmt(myfile)
    cf <- data.frame(geneset = d$geneset.names,
                     num_genes=NA,
                     genes=NA,
                     genes_with_alias=NA)
    for(i in 1:length(d$genesets)){
      cat(i,"\n")
      cf$num_genes[i] <- length(d$genesets[[i]])
      cf$genes[i] <- paste0(d$genesets[[i]],collapse = ",")
      myhgnc <- unique(genes[genes$geneSymbol%in%d$genesets[[i]],]$HGNC)
      cf$genes_with_alias[i] <- paste0(unique(genes[genes$HGNC%in%myhgnc,]$geneSymbol),collapse = ",")
      rm(myhgnc)
    }
    out <- rbind(out,cf)
    rm(d,cf)
  }
  
  write.table(out[1:50,],file = "msigdb_signatures.txt",quote = F,sep = "\t",row.names = F)
  
  
}

##--- general tests
if(F){
  source("utilities.R",local = TRUE,encoding = "UTF-8")
  source("vars.R",local = TRUE,encoding = "UTF-8")
  
  ##-- from rds
  so <- local({
    sa <- readRDS("/home/rstudio/data/SCS/TM_BoneMarrow_SmartSeq.rds")
    si_assays<- Seurat::Assays(sa)
    so <- CreateSeuratObject(counts = sa@assays[[1]]@data,
                             assay = si_assays[1],
                             meta.data = sa@meta.data)
    if(length(si_assays)>1){
      message(' ..adding additional assay information to the Seurat object\n')
      for(i in 2:length(si_assays)){
        so[[si_assays[i]]] <- CreateAssayObject(counts = sa@assays[[i]]@data)
      }
    }
    scvis <- sa@reductions
    if(length(scvis)>1) message('Multiple data reductions found. Will populate them.\n')
    for(i in 1:length(scvis)){
      message(" ..adding coordinates: ",gsub("^X_","",names(scvis)[i]),"\n")
      so[[gsub("^X_","",names(scvis)[i])]] = CreateDimReducObject(embeddings = scvis[[i]]@cell.embeddings, key = gsub("^X_","",names(scvis)[i]))
    }
    so
  })
  so@active.ident <- factor(so$free_annotation)
  
  ##-- from h5ad file
  so <- local({
    ad <- anndata::read_h5ad(filename = "/home/rstudio/data/SCS/HS_healthyCholangiocytes_10x.h5ad")
    meta <- ad$obs
    meta$V1 <- meta$V2 <- NULL
    mat <- ad$X
    qq <- dimnames(mat)
    mat <- as(mat,'matrix.csr')
    mat <- as(mat,'dgCMatrix')
    dimnames(mat) <- qq
    rm(qq)
    mat <- t(mat)
    cat(dim(mat)[1],"\n")
    so = CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0,meta.data = meta)
    scvis <- ad$obsm
    if(length(scvis)>1) message('Multiple data reductions found. Will populate them.\n')
    for(i in 1:length(scvis)){
      message(" ..adding coordinates: ",gsub("^X_","",names(scvis)[i]),"\n")
      so[[gsub("^X_","",names(scvis)[i])]] = CreateDimReducObject(embeddings = scvis[[i]], key = gsub("^X_","",names(scvis)[i]))
    }
    so
  })
  
  
  # subsample
  sampled.cells <- sample(x = colnames(so), 
                          size = 100, replace = F)
  so <- subset(so, cells = sampled.cells)
  rm(sampled.cells)
  
  ##-- meta
  meta <- local({
    coords <- data.frame(SAMPID=rownames(so@meta.data))
    for(f in names(so@reductions)){
      qq <- as.data.frame(so@reductions[[f]]@cell.embeddings)
      qq <- qq[,1:2]
      names(qq) <- c('V1','V2')
      names(qq) <- paste0(names(qq),'_',f)
      coords <- cbind(coords, qq)
      rm(qq)
    }
    my.meta <- so@meta.data
    my.meta <- cbind(my.meta,coords)
    my.meta
  })
  
  
  study <- data.frame(Database="abc",
                      Description=NA,
                      SampleSize=NA,
                      TISSUES=NA,
                      DISEASE_Variable=NA, 
                      PMID=NA,
                      GEO=NA,
                      STATUS='Public',
                      RATING='High',
                      CONTVAR=NA,
                      CATVAR='free_annotation')
  
  seurat2sqlite(so=so,
                si_study=study,
                si_reduction='umap',
                si_compute_cellmarker = F,
                si_cell_markers=NULL,
                si_cell_col = NULL,
                si_compute_diseasemarker = F,
                si_disease_markers=NULL,
                si_disease_col=NULL,
                si_celltype='free_annotation',
                si_donor = 'mouse.id',
                signature_file = NULL,
                db_address='asd',
                OUTDIR=NULL)
  
  
  si_celltype='free_annotation'
  si_donor = 'mouse.id'
  my_assay='RNA'
  signature_file <- read.delim("data/msigdb_signatures.txt",header=F,stringsAsFactors = F)
  si_reduction='umap'
  si_study=study
  si_compute_cellmarker = F
  si_cell_markers=NULL
  si_cell_col = NULL
  si_compute_diseasemarker = F
  si_disease_markers=NULL
  si_disease_col=NULL
  db_address='asd'
  OUTDIR=NULL

  
  ## quick checks
  DefaultAssay(so) = 'RNA'
  group.by = c(si_celltype,si_donor)
  
  m <- so@assays$RNA@data
  n <- m[,colnames(m)[1:20]]
  a <- apply(n,1,mean)
  
  m <- FetchData(so,vars = c(si_celltype,si_donor))
  names(m) <- c('cell_type','donor')
  m$id <- paste0(m$cell_type,"|",m$donor)
  split(m,m$id)

  
  m <- m[which(rowSums(is.na(m)) == 0), , drop = F]
  for (i in 1:ncol(m)) {
    m[, i] <- as.factor(m[, i])
  }
  num.levels <- sapply(
    X = 1:ncol(m),
    FUN = function(i) {
      length(levels(m[, i]))
    }
  )
  
  apply(so@assays$RNA@data,1,sd)
  
  
  category.matrix <- sparse.model.matrix(object = as.formula(
    object = paste0(
      '~0+',paste0("m[,",1:length(group.by),"]",collapse = ":"
      )
    )))
  colsums <- colSums(x = category.matrix)
  category.matrix <- category.matrix[, colsums > 0]
  colsums <- colsums[colsums > 0]
  colnames(category.matrix) <- sapply(
    X = colnames(x = category.matrix),
    FUN = function(name) {
      name <- gsub(pattern = "m\\[, [1-9]*\\]", replacement = "", x = name)
      return(paste0(rev(x = unlist(x = strsplit(x = name, split = ":"))), collapse = "_"))
    })
  
  
  m <- Seurat::AverageExpression(object = so,assays = my_assay,
                            group.by = c(si_celltype,si_donor))[[1]]
  
}

##--- test output
if(F){
  db_address = "/home/rstudio/data/SCS/scRNA_RNA_Hs_TabulaSapiensStroma_10x.db"
  connSc <- RSQLite::dbConnect(RSQLite::SQLite(),paste0(db_address))
  RSQLite::dbListTables(connSc)
  RSQLite::dbDisconnect(connSc)
  
  study = "RNA_Hs_TabulaSapiensImmune_10x"
  sc_study <- RSQLite::dbGetQuery(conn = connSc,paste0('SELECT * FROM ',study,"_study"))
  m <- RSQLite::dbGetQuery(conn = connSc,paste0('SELECT * FROM ',study,"_metaFeatures"))
  n <- RSQLite::dbGetQuery(conn = connSc,paste0('SELECT * FROM ',study,"_FeatureSummary"))
  
  
  #c("author_cell_type","donor_uuid","Phase","cell_type","disease",
  #  "development_stage","sample_preservation_method","sex","author_cluster")
  Categorical_Vars <- "author_cell_type,donor_uuid,Phase,cell_type,disease,development_stage,sample_preservation_method,sex,author_cluster"
  Celltype <- "author_cell_type"
  qq <- unlist(strsplit(as.character(Categorical_Vars),","))
  qq <- ifelse(qq==Celltype,as.character("SVcell_type"),as.character(qq))
  if(sum(qq%in%names(so@meta.data))!=length(qq) ){
    #cat(paste0(qq,collapse = ","),"\n",
    #    paste0(names(so@meta.data),","),"\n")
    stop(" could not find provided categorical variables in the metadata file\n") 
  }
  if(!'SVcell_type'%in%qq){
    stop(" could not find column with name 'cell_type' in the metadata file\n") 
  }
  Categorical_Vars = qq
  message('Categorical variables in the metadata: ',Categorical_Vars,"\n\n")
  
   
  
  study = "RNA_test"
  db_address = "/home/rstudio/data/SCS/scRNA_RNA_test.db"
  connSc <- RSQLite::dbConnect(RSQLite::SQLite(),paste0(db_address))
  RSQLite::dbListTables(connSc)
  
  sc_study <- RSQLite::dbGetQuery(conn = connSc,paste0('SELECT * FROM ',study,"_study"))
  m <- RSQLite::dbGetQuery(conn = connSc,paste0('SELECT * FROM ',study,"_metaFeatures"))
  n <- RSQLite::dbGetQuery(conn = connSc,paste0('SELECT * FROM ',study,"_FeatureSummary"))
  o <- RSQLite::dbGetQuery(conn = connSc,paste0('SELECT * FROM ',study,"_FeatureSummaryKeys"))
  
  RSQLite::dbDisconnect(connSc)
  
}

##---------------------------------------------------
## single cell studies over time
##---------------------------------------------------
if(F){
  cf <- data.frame(year=2003:2021,
                   publications=c(0,1,1,1,1,6,30,123,357,727,1186,1617,3204,4919,6752,9253,11824,14623,16399))
  
  ggplot(cf,aes(x=year,y=publications,color=publications,size=(publications)))+
    geom_text(aes(label=publications))+
    theme_bw()+
    scale_color_gradient(low = '#fc4e2a',high = '#800026')+
    theme(legend.position = 'none',
          axis.text = element_text(color='black',size=14),
          axis.title = element_text(size=14))
  
  
}


##---------------------------------------------------
## CCF single cell: Bad data! 
##---------------------------------------------------
if(F){
  rm(list=ls())
  gc()
  setwd("~/scripts/SciViewerDev/SciViewIn/")
  source("/home/rstudio/scripts/SciViewerDev/SciViewIn/utilities.R",local = TRUE,encoding = "UTF-8")
  source("/home/rstudio/scripts/SciViewerDev/SciViewIn/vars.R",local = TRUE,encoding = "UTF-8")
  setwd("~/data/SCS/")
  
  filepath = "/home/rstudio/data/SCS/HS_CrohnDiseaseCCF_1_10x.h5ad"
  
  ad <- anndata::read_h5ad(filename = filepath)
  meta <- ad$obs
  meta$V1 <- meta$V2 <- NULL
  mat <- ad$raw$X
  qq <- dimnames(ad$X)
  mat <- as(mat,'matrix.csr')
  mat <- as(mat,'dgCMatrix')
  dimnames(mat) <- qq
  rm(qq)
  mat <- t(mat)
  cat(dim(mat)[1],"\n")
  so = CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0,meta.data = meta)
  scvis <- ad$obsm
  for(i in 1:length(scvis)){
    message(" ..adding coordinates: ",gsub("^X_","",names(scvis)[i]),"\n")
    so[[gsub("^X_","",names(scvis)[i])]] = CreateDimReducObject(embeddings = scvis[[i]], key = gsub("^X_","",names(scvis)[i]))
  }
  meta = so@meta.data
  meta$condition <- gsub("CD2586|CD2587","",meta$patient)
  meta$condition <- ifelse(meta$condition=="i","inflammed",meta$condition)
  meta$condition <- ifelse(meta$condition=="ni","non-inflammed",meta$condition)
  meta$condition <- ifelse(meta$condition=="s","stricture",meta$condition)
  meta$patient <- gsub("i|s|ni","",meta$patient)
  so@meta.data = meta
  saveRDS(so,file="HS_CrohnDiseaseCCF_1_10x.rds")
  
  so <- readRDS(file="HS_CrohnDiseaseCCF_1_10x.rds")
  meta <- so@meta.data
  names(meta)[7] <- "layer"
  meta$condition_layer <- paste0(meta$condition,"_",meta$layer)
  meta$condition_layer_v2 <- paste0(meta$condition,"_",meta$layer)
  meta$condition_layer_v2 <- gsub("_LP1|_LP2","_LP",meta$condition_layer_v2)
  table(meta$condition_layer_v2)
  table(so$cell_type)
  so@meta.data = meta
  saveRDS(so,file="HS_CrohnDiseaseCCF_1_10x.rds")
  
  
}



##---------------------------------------------------
## batch run
##---------------------------------------------------
if(F){
  rm(list=ls())
  gc()
  setwd("~/scripts/SciViewerDev/SciViewIn/")
  source("/home/rstudio/scripts/SciViewerDev/SciViewIn/utilities.R",local = TRUE,encoding = "UTF-8")
  source("/home/rstudio/scripts/SciViewerDev/SciViewIn/vars.R",local = TRUE,encoding = "UTF-8")
  setwd("~/data/SCS/")
  
  
  files <- Sys.glob("*_params.txt")
  files <- files[c(5,14,7)]
  for(f in files){
    cf = read.delim(file = f,header = T,stringsAsFactors = F,as.is = T)
    cat(paste0(cf[cf$field=="inptfile_path",]$value,"\n"))
    cf$value <- as.character(cf$value)
    cf$field <- as.character(cf$field)
    so <- readInpscRNAseq(filepath=cf[cf$field=="inptfile_path",]$value,intractv=F)
    tryCatch({
      write.scstudy2.sqlitedb(so = so,
                              db_address = cf[cf$field=="db",]$value,
                              StudyName = cf[cf$field=="studyname",]$value,
                              Celltype = cf[cf$field=="cell_type",]$value,
                              Reduction_map = cf[cf$field=="redmap",]$value,
                              Donors_VariableName = cf[cf$field=="donor",]$value,
                              Continuous_Vars = cf[cf$field=="contvars",]$value,
                              Categorical_Vars = cf[cf$field=="catvars",]$value,
                              Marker_Calc = cf[cf$field=="simarker",]$value,
                              Marker_Covariates = cf[cf$field=="covariates",]$value,
                              Marker_Precomputed = cf[cf$field=="si_markerfile_path",]$value,
                              DE_Calc = cf[cf$field=="de",]$value,
                              Disease_VariableName = cf[cf$field=="disease",]$value,
                              DE_Precomputed = cf[cf$field=="si_defile_path",]$value,
                              StudyDescr = cf[cf$field=="descr",]$value,
                              Tissue = cf[cf$field=="tissue",]$value,
                              Organism = cf[cf$field=="svorg",]$value,
                              Disease = cf[cf$field=="svdisease",]$value,
                              ShortDescr = cf[cf$field=="shotdescr",]$value,
                              PMID = cf[cf$field=="pmid",]$value,
                              GEO = cf[cf$field=="geo",]$value,
                              StudyStatus = cf[cf$field=="status",]$value,
                              StudyRating = cf[cf$field=="rating",]$value,  
                              OUTDIR = paste0(dirname(cf[cf$field=="inptfile_path",]$value),"/final"),
                              signature_file=NULL,
                              intractv=F)
    },error=function(e){ return(NA)})
    rm(cf,so,f)
    gc()
    
  }
  
}










