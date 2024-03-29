source("vars.R",local = TRUE,encoding = "UTF-8")

##------------------------------------------------------------------------
## Input functions
##------------------------------------------------------------------------

remove_dupcolumns <- function(temp){
  temp[, !duplicated(colnames(temp))]
}


matrix.reformat <- function(mat,meta,intractv=F){
  
  require(data.table)|| stop("Can't find data.table package. Please install and restart! \n")
  require(Matrix)|| stop("Can't find Matrix package. Please install and restart! \n")
  
  if(intractv==T){
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = "reformatting matrix", value = 0)
  }
  if(!'SAMPID'%in%names(meta)){
    if(intractv==T){
      scProg$set(message = "SAMPID column is not found in the metadata file. Please make sure to label the cell id column as SAMPID", value = 0.1)
    }
    stop('SAMPID column is not found in the metadata file. Please make sure to label the cell id column as SAMPID\n')
  }
  if(length(which(is(mat)=="sparseMatrix"))==0){
    if(intractv==T){
      scProg$set(message = "provided input matrix is not in sparsematrix format!\n Converting data to sparse matrix format.", value = 0.2)
    }
    warning(' provided input matrix is not in sparsematrix format!\n Converting data to sparse matrix format.')
    mat <- as.matrix(mat)
    mat <- as(mat,'dgCMatrix')
  }
  genes <- data.frame(gene_id=rownames(mat),geneSymbol=rownames(mat))
  
  
  meta <- meta[meta$SAMPID%in%colnames(mat),]
  mat <- mat[,colnames(mat)%in%meta$SAMPID]
  if(ncol(mat)!=sum(colnames(mat)%in%meta$SAMPID)){
    if(intractv==T){
      scProg$set(message = "Provided number of cells in the count matrix do not match to those in the metadata file", value = 0.3)
    }
    stop('Provided number of cells in the count matrix do not match to those in the metadata file')
  }
  
  if(intractv==T){
    scProg$set(message = paste0('Matrix dimensions: ',paste0(dim(mat),collapse = ",")), value = 0.4)
    scProg$set(message = paste0('Generating indexes'), value = 0.5)
  }
  message('Matrix dimensions: ',paste0(dim(mat),collapse = ","),"\n")
  message('Generating indexes..\n')
  rnames <- rownames(mat)
  rnames <- data.frame(geneID=rnames,id=1:length(rnames))
  cnames <- colnames(mat)
  cnames <- data.frame(cellID=cnames,id=1:length(cnames))
  mat <- as.data.frame(summary(mat))
  
  if(intractv==T){
    scProg$set(message = paste0('Regrouping the data (this step takes relatively long time to finish, please wait)'), value = 0.7)
  }
  message("Regrouping the data..\n")
  message('  -->this step takes relatively long time to finish, please wait..')
  mat <- data.table(mat)
  sumfun <- function(x){ paste0(x,collapse = ",") }
  mat <- mat[, lapply(.SD, sumfun), by=list(i) ]
  mat <- as.data.frame(mat)
  names(mat)[1:3] <- c("row_index","col_index",'value')
  
  ## match the gene-dataframe with matrix index
  rnames <- rnames[match(mat$row_index,rnames$id),]
  genes <- genes[match(rnames$geneID,genes$gene_id),]
  if(nrow(mat)!=sum(genes$gene_id==rnames$geneID)){
    if(intractv==T){
      scProg$set(message = paste0('error, the rows do not match!'), value = 0.9)
    }
    stop('error, the rows do not match! \n')
  }
  
  ## binding the matrix information
  mat <- cbind(genes,mat)
  rownames(mat) <- NULL
  
  ## re-arrange meta dataframe rows to match original matrix index
  meta <- meta[match(cnames$cellID,meta$SAMPID),]
  meta$col_index <- cnames$id
  
  if(intractv==T){
    scProg$set(message = "Done!", value = 1)
  }
  return(list(A=mat,B=meta))
  
}


seurat2expression.summary <- function(so,si_celltype=NA,
                                      si_donor=NA,my_assay='RNA',intractv=F){
  
  if(intractv==T){
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = "computing expression summaries..", value = 0)
  }
  if(is.na(si_celltype)){
    if(intractv==T){
      scProg$set(message = "No celltype information is specified for expression aggregation", value = 0.1)
    }
    stop('No celltype information is specified for expression aggregation\n')
  }
  if(intractv==T){
    scProg$set(message = "by default the function averages the gene expression as normalized by the user", value = 0.2)
  }
  message("  by default the function averages the gene expression as normalized by the user\n")
  
  if(is.na(si_donor)){
    x <- AverageExpression(object = so,features = rownames(so),return.seurat = FALSE,
                           group.by = c(si_celltype),slot = "data",verbose = TRUE)
    }else{
      x <- AverageExpression(object = so,features = rownames(so),return.seurat = FALSE,
                             group.by = c(si_celltype,si_donor),slot = "data",verbose = TRUE)
  }
  if(intractv==T){
    scProg$set(message = paste0('  ..computed average expression values'), value = 0.8)
  }
  message('  ..computed average expression values')
  
  x <- x[[1]] 
  x <- as.data.frame(x)
  if(is.na(si_donor)){
    names(x) <- paste0(names(x),'_A') 
    mycomb <- data.frame(Var1=levels(as.factor(so@meta.data[,si_celltype])),Var2='A')
  }else{
    mycomb <- expand.grid(levels(as.factor(so@meta.data[,si_celltype])),
                          levels(as.factor(so@meta.data[,si_donor])),stringsAsFactors = F)
  }
  mycomb$id <- paste(mycomb$Var1,mycomb$Var2,sep = "_")
  mycomb <- mycomb[match(names(x),mycomb$id),]
  if(sum(names(x)==mycomb$id)!=ncol(x)) stop('  ..the cell/donor column combination is not fully captured while generating expression summary')
  mycomb$code <- paste0("ab",1:nrow(mycomb))
  names(x) <- mycomb$code
  x <- cbind(feature=rownames(x),x)
  rownames(x) = NULL
  names(mycomb)[1:2] <- c('cell_type','donor')
  
  #expn.summaryperfeature <- c()
  #m <- FetchData(so,vars = c(si_celltype,si_donor))
  #if(is.null(si_donor)) m$donor <- "A"
  #names(m) <- c('cell_type','donor')
  #m$id <- paste0(m$cell_type,"|",m$donor)
  #m <- split(m,m$id)
  #my.data <- so@assays[[my_assay]]@data
  #for(i in 1:length(m)){
  #  if(intractv==T){
  #    scProg$set(message = paste0('  ..working on expression group: ',i,"/",length(m),' (',names(m)[i],')'), value = 0)
  #  }
  #  message('  ..working on expression group: ',i,"/",length(m),' (',names(m)[i],')')
  #  x <- my.data[,rownames(m[[i]])]
  #  if(nrow(m[[i]])>1){
  #    y <- data.frame(
  #      cell_type=unique(m[[i]]$cell_type),
  #      donor = unique(m[[i]]$donor),
  #      feature = rownames(x),
  #      norm_avg_priorLT= apply(x,1,function(x) { log2(mean(2^x,na.rm=T)) }),
  #      norm_avg= apply(x,1,mean,na.rm=T),
  #      norm_sd = apply(x,1,sd,na.rm=T),
  #      norm_gmean = 0
  #    )
  #    expn.summaryperfeature <- rbind(expn.summaryperfeature,y)
  #  }else{
  #    y <- data.frame(
  #      cell_type=unique(m[[i]]$cell_type),
  #      donor = unique(m[[i]]$donor),
  #      feature = names(x),
  #      norm_avg_priorLT= x,
  #      norm_avg= x,
  #      norm_sd = 0,
  #      norm_gmean = 0
  #    )
  #    expn.summaryperfeature <- rbind(expn.summaryperfeature,y)
  #  }
  #  rm(x,y)
  #}
  #rm(my.data,m,i)
  #gc()
  #expn.summaryperfeature <- expn.summaryperfeature[abs(expn.summaryperfeature$norm_avg)>0,]
  
  if(intractv==T){
    scProg$set(message = "Done!", value = 1)
  }
  return(  list(expn.summaryperfeature=x,mykeys=mycomb) )
}


seurat2.signatures <- function(sa,signature_file=NULL,my_assay='RNA',intractv=F){
  if(is.null(signature_file)){
    message(" ..gene signature file does not exits. Skipping signature calculations.\n")
    return(NULL)
  }
  
  DefaultAssay(sa) = my_assay
  signature_file$total <- signature_file$captured <- 
    lapply(1:nrow(signature_file),function(i){
    x <- unlist(strsplit(signature_file$V3[i],","))
    y <- sum(x%in%rownames(sa))
    return(paste0(y,"|",length(x)))
  }) %>% unlist()
  signature_file$total <- gsub("\\d*\\|","",signature_file$total) %>% as.numeric()
  signature_file$captured <- gsub("\\|\\d*","",signature_file$captured) %>% as.numeric()
  signature_file$prop <- signature_file$captured/signature_file$total
  signature_file <- signature_file[signature_file$prop>0.8 & signature_file$captured>3,]
  if(isTRUE(nrow(signature_file)==0)){
    message(" ..no gene signatures were computed on the data.\n")
    return(NULL)
  }
  message(' calculating the signatures for ', nrow(signature_file),' gene sets!')
  
  l <- lapply(1:nrow(signature_file),function(i){
     unlist(strsplit(signature_file$V3[i],","))
  })
  names(l) <- signature_file$V1  
  meta = sa@meta.data
  meta = meta[,1:2]
  sa@meta.data=meta
  sa <- AddModuleScore(sa,features = l,name=names(l))
  meta=sa@meta.data
  meta= meta[,3:ncol(meta)]
  colnames(meta) = signature_file$V1
  
  #m  <- data.frame(SAMPID=rownames(meta),
  #                 signature.value = apply(meta,1,function(x) paste0(x,collapse = ",")))
  return(list(sig=meta,sigfile=signature_file))
}


writedb <- function(tableName=NULL,datafrm=NULL,db_address,OUTDIR=NULL,overwrite=T){
  if(!is.null(tableName) & !is.null(datafrm)){
    if(!is.null(OUTDIR)){
      connSc <- RSQLite::dbConnect(RSQLite::SQLite(),paste0(OUTDIR,"/",db_address))
    }else{
      connSc <- RSQLite::dbConnect(RSQLite::SQLite(),paste0(db_address))
    }
    datafrm <- remove_dupcolumns(datafrm)
    RSQLite::dbWriteTable(connSc,tableName,datafrm,overwrite=overwrite)
    RSQLite::dbListTables(connSc)
    RSQLite::dbDisconnect(connSc)
  }
}


seurat2sqlite <- function(so=NULL,si_study=NULL,
                          si_reduction=NULL,
                          si_compute_cellmarker = F,si_cell_markers=NULL,si_cell_col = NULL,
                          si_compute_diseasemarker = F,si_disease_markers=NULL,si_disease_col=NULL,si_celltype=NULL,
                          si_donor = NULL,
                          db_address=NULL,OUTDIR=NULL,
                          signature_file=NULL,intractv=F){
  
  if(intractv==T){
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = "assessing dependancies..", value = 0)
  }
  
  require(Seurat) || stop("Can't find Seurat package. Please install and restart! \n")
  require(RSQLite) || stop("Can't find RSQLite package. Please install and restart! \n")
  #-- dependencies
  if(is.null(so)) stop(" Please provide Seurat object for recording the study\n")
  if(is.null(si_study)) stop(" Please provide study name for recording the study\n")
  if(is.null(db_address)) stop(" Please provide DB name for recording the study\n")
  
  
  #-- collect variables
  si_assays<- Seurat::Assays(so)
  si_coords <- names(so@reductions)
  if(!is.null(si_reduction) & si_reduction %in% names(so@reductions)){
    si_default_coord <- si_reduction
  }else{
    si_default_coord <- ifelse(grep('^umap$|^UMAP$',si_coords),
                               si_coords[grep('^umap$|^UMAP$',si_coords)],
                               si_coords[1])
  }
  si_study$COORDS <- paste0(si_coords,collapse = ",")
  si_study$DEFAULT_COORDS <- si_default_coord
  si_study$cell_type <- si_celltype 
  if(is.null(si_donor)){
    so$tmpdonor <- 1
    si_donor = "tmpdonor"
  } 
  si_study$donor <- si_donor
  si_study$FeatureCount <- dim(so)[1] 
  si_study$CellCount <- dim(so)[2] 
  
  #-- data normalization etc.
  if(length(si_assays)>1){
    warning('The provided Seurat object contains multiple assay objects.Please make sure that for non-RNA assays you have appropriately normalized the data.RNA assay object will be tested by default and log-normalized in case non-normalized data is found.\n')
  }
  message(paste0("IMPORTANT: Range of the provided RNA data: ",range(so[[si_assays[1]]]@data),"\n"))
  if(range(so[[si_assays[1]]]@data)[2]>25){
    if(intractv==T){
      scProg$set(message = "Normalizing the RNA data..", value = 0.3)
    }
    message('Normalizing the RNA data\n')
    so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)  
  }
  
  #-- add additional qc estimates
  #tryCatch({
  #  if(!'percent.mt'%in%names(so@meta.data)) so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  #  if(!'percent.RPS'%in%names(so@meta.data)) so[["percent.RPS"]] <- PercentageFeatureSet(so, pattern = "^RPS")
  #  if(!'percent.RPL'%in%names(so@meta.data)) so[["percent.RPL"]] <- PercentageFeatureSet(so, pattern = "^RPL")
  # },error = function(e){
  #   scProg$set(message = "could not add additional QC parameters..", value = 0.3)
  # })
  
  #-- create metadata file
  my.meta <- local({
    coords <- data.frame(SAMPID=rownames(so@meta.data))
    qq <- as.data.frame(so@reductions[[si_default_coord]]@cell.embeddings)
    qq <- qq[,1:2]
    names(qq) <- c('V1','V2')
    coords <- cbind(coords, qq)
    rm(qq)
    for(f in si_coords){
      qq <- as.data.frame(so@reductions[[f]]@cell.embeddings)
      qq <- qq[,1:2]
      names(qq) <- c('V1','V2')
      names(qq) <- paste0(names(qq),'_',f)
      coords <- cbind(coords, qq)
      rm(qq)
    }
    my.meta <- so@meta.data
    my.meta <- cbind(my.meta,coords)
    my.meta <- remove_dupcolumns(my.meta)
    my.meta
  })
  if(intractv==T){
    scProg$set(message = "created metadata..", value = 0.3)
  }
  
  #-- loop over assays 
  for(my_assay in si_assays){
    message("working on assay:", my_assay,"\n")
    DefaultAssay(so) <- my_assay
    si_cell_markers_tmp <- c()
    si_disease_markers_tmp <- c()
    expn.summaryperfeature <- c()
    signature.list <- c()
    si_study$assay <- my_assay
    
    #--------------------------------
    #--(1) study dataframe
    #--------------------------------
    if(isTRUE(nrow(si_study)>0)){
      if(intractv==T){
        scProg$set(message = "Writing down study metadata", value = 0.35)
      }
      message('Writing down study metadata\n')
      StudyName = paste0(my_assay,"_",si_study$Database)
      si_study$Database = StudyName
      if(length(grep("^scRNA_",db_address))==0){
        db_address <- paste0("scRNA_",my_assay,"_",db_address)
      }
      writedb(tableName=paste0(StudyName,"_study"),datafrm=si_study,db_address=db_address,OUTDIR=OUTDIR,overwrite=T)
    }else(
      stop(" single cell study dataframe is not found. Exiting!")
    )
    
    #--------------------------------
    #--(2) matrix reformatting
    #--------------------------------
    if(intractv==T){
      scProg$set(message = "Reformatting the matrix for writing down", value = 0.4)
    }
    message('Reformatting the matrix for writing down\n')
    l<- matrix.reformat(mat=so@assays[[my_assay]]@data,
                        meta=my.meta,intractv = intractv)
    if(isTRUE(nrow(l[[1]])>0) & isTRUE(nrow(l[[2]])>0)){
      writedb(tableName=paste0(StudyName,"_metaFeatures"),datafrm=l[[2]],db_address=db_address,OUTDIR=OUTDIR,overwrite=T)
      writedb(tableName=paste0(StudyName,"_data"),datafrm=l[[1]],db_address=db_address,OUTDIR=OUTDIR,overwrite=T)
    }
    rm(l)
    gc()
    
    
    #--------------------------------
    #--(3) compute gene expression summary
    #--------------------------------
    message(' computing the expression summary for each feature. \n')
    if(intractv==T){
      scProg$set(message = "computing the expression summary for each feature...", value = 0.5)
    }
    expn.summaryperfeature <- seurat2expression.summary(so=so,
                                                        si_celltype = si_celltype,
                                                        si_donor = si_donor,
                                                        my_assay = my_assay,
                                                        intractv=intractv)
    if(isTRUE(nrow(expn.summaryperfeature[[1]])>0)){
      writedb(tableName=paste0(StudyName,"_FeatureSummary"),datafrm=expn.summaryperfeature[[1]],db_address=db_address,OUTDIR=OUTDIR,overwrite=T)
      writedb(tableName=paste0(StudyName,"_FeatureSummaryKeys"),datafrm=expn.summaryperfeature[[2]],db_address=db_address,OUTDIR=OUTDIR,overwrite=T)
    }
    rm(expn.summaryperfeature)
    gc()
    
    
    #--------------------------------
    #--(4) signature calculations
    #--------------------------------
    if(!is.null(signature_file) & my_assay=='RNA'){
      message(' Calculating gene signatures using aggregation method.\n')
      if(intractv==T){
        scProg$set(message = "Calculating gene signatures using aggregation method....", value = 0.6)
      }
      signature.list <- seurat2.signatures(sa = so,
                              signature_file = signature_file,
                              my_assay = my_assay)
    }
    if(isTRUE(length(signature.list)>0)){
      writedb(tableName=paste0(StudyName,"_SignatureMatrix"),datafrm=signature.list[[1]],db_address=db_address,OUTDIR=OUTDIR,overwrite=T)
      writedb(tableName=paste0(StudyName,"_SignatureFile"),datafrm=signature.list[[2]],db_address=db_address,OUTDIR=OUTDIR,overwrite=T)
    }
    rm(signature.list)
    gc()
    
    #--------------------------------
    #--(5) compute disease markers
    #--------------------------------
    if(si_compute_diseasemarker==T &  !is.null(si_disease_col) & !is.null(si_celltype)){
      my.disease <- levels(as.factor(so@meta.data[,si_disease_col]))
      my.celltypes <- levels(as.factor(so@meta.data[,si_celltype]))
      if(length(my.disease)<2) stop(" the disease/exp group doesn't have contrast\n")
      if(length(si_celltype)>1) stop(" please define unique cell type annotation column\n")
      if(length(si_disease_col)>1) stop(" please define unique disease/experimental condition annotation column\n")
      if(!si_disease_col%in%names(so@meta.data)) stop(' the disease/experimental condition column is not found in the metadata\n')
      my.tests <- expand.grid(my.disease,my.disease,stringsAsFactors = F)
      my.tests <- my.tests[my.tests$Var1>my.tests$Var2,]
      message(" Performing the differential expression tests using provided disease/experiment contrast.\n")
      message("  ..total ",nrow(my.tests)," contrasts are found. Starting the execution!\n")
      message("  ..this can take a while.\n")
      if(intractv==T){
        scProg$set(message = "  ..computing DEG", value = 0.8)
      }
      so$temp_celltype <- so@meta.data[,si_celltype]
      Idents(so) <- so@meta.data[[si_disease_col]]
      for(i in 1:nrow(my.tests)){
        sa <- subset(x = so, idents = c(my.tests$Var1[i],my.tests$Var2[i]))
        for(f in my.celltypes){
          sb <- tryCatch({
            subset(x = sa, subset = temp_celltype == f)
          },error = function(e){
            return(c())
          })
          if(is.null(sb)) next
          if(intractv==T){
            scProg$set(message = paste0('  ..working on celltype: ',f ,"  (",my.tests$Var1[i],"-", my.tests$Var2[i],")\n"), value = 0.9)
          }
          message(paste0('  ..working on celltype: ',f ,"  (",my.tests$Var1[i],"-", my.tests$Var2[i],")\n"))
          q <- tryCatch({
            plan("multicore", workers = MCORES)
            FindMarkers(sb, ident.1 = my.tests$Var1[i], 
                        ident.2 = my.tests$Var2[i],
                        print.bar = T,
                        min.pct = 0.25, 
                        logfc.threshold = 0.25,
                        return.thresh = 0.1,
                        max.cells.per.ident = 2000,
                        min.cells.group = 3,
                        test.use = 'wilcox')
          },error=function(e){
            return(c())
          })
          if(!is.null(q)){
            names(q) <- c("P.Value", "logFC", "pct.1", "pct.2", 
                          "adj.P.Val")
            q$Test <- paste0(my.tests$Var1[i],"-", my.tests$Var2[i])
            q$Tag <- f
            q$geneSymbol <- rownames(q)
            q$Method <- 'Wilcox'
            q$Assay <- my_assay
            q$AveExpr <- NA
            q$t <- NA
            q$B <- NA
            q <- q[,c("logFC", "AveExpr", "t","P.Value", "adj.P.Val",
                      "B", "Test","Tag","geneSymbol",'Method','Assay',
                      "pct.1","pct.2")]
            si_disease_markers_tmp <- rbind(si_disease_markers_tmp,q)
          }
          rm(q,sb)
        }
      }
      so$temp_celltype <- NULL
      
    }else if(si_compute_diseasemarker==F & !is.null(si_disease_col) & !is.null(si_celltype)){
      mycols <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Test","Tag","geneSymbol")
      cf <- NULL
      ext <- tools::file_ext(si_disease_markers)
      if(ext=='txt'){
        cf = read.delim(si_disease_markers,as.is = T,colClasses = "#",sep="\t")
      }else if(ext == 'csv'){
        cf = read.csv2(si_disease_markers,as.is = T,colClasses = "#")
      }
      if(sum(mycols%in%names(cf))!=length(mycols)){
        warning("One or more columns are missing in the provided disease marker file. Please ensure that correct file is provided. Skipping this file and moving on!")
      }else{
        mycols <- append(mycols,
                         names(cf)[!names(cf)%in%mycols]
        )
        si_disease_markers_tmp <- cf[,mycols]
      }
      rm(cf)
    }
    if(isTRUE(nrow(si_disease_markers_tmp)>0)){
      writedb(tableName=paste0(StudyName,"_DEG"),datafrm=si_disease_markers_tmp,db_address=db_address,OUTDIR=OUTDIR,overwrite=T)
    }
    rm(si_disease_markers_tmp)
    gc()
    plan("sequential")
    
    
    #--------------------------------
    #--(6) compute cell type markers
    #--------------------------------
    if(intractv==T){
      scProg$set(message = paste0(" Computing markers for the variable/s:", si_cell_col,"\n"), value = 0.7)
    }
    si_cell_col = unlist(strsplit(as.character(si_cell_col),","))
    if(si_compute_cellmarker==T & !is.null(si_cell_col)){
      if(is.null(si_cell_col)) stop(" Column name/s of the variable against which markers are to be computed is not provided!\n")
      if(sum(si_cell_col%in%names(so@meta.data))!=length(si_cell_col)){
        cat(paste0(si_cell_col,collapse = ","),"\n")
        cat(paste0(names(so@meta.data),collapse = ","),"\n")
        cat(sum(si_cell_col%in%names(so@meta.data)),"  ",length(si_cell_col),"\n")
        stop(' one or more marker column names are not present in the Seurat object\n')
      }
      message(" Computing markers for the variable/s:", si_cell_col,"\n")
      if(length(si_cell_col)>1){
        message(' -> Note: multiple varialbes are provided for computing the markers. Will proceed!\n')
      }
      for(mycolm in si_cell_col){
        message("  ..computing markers for the variable:", mycolm,"\n")
        if(intractv==T){
          scProg$set(message = paste0("  ..computing markers for the variable:", mycolm,"\n"), value = 0.75)
        }
        Idents(so) <- so@meta.data[[mycolm]]
        plan("multicore", workers = MCORES)
        my.markers <- tryCatch({FindAllMarkers(so,  
                                               min.pct = 0.25, 
                                               logfc.threshold = 0.25,
                                               return.thresh = 0.1,
                                               min.diff.pct = 0.25,
                                               max.cells.per.ident = 500,
                                               test.use = 'wilcox')
        },error=function(e){
          return(c())
        })
        if(!is.null(my.markers)){
          if(ncol(my.markers)==7){
            names(my.markers) <- c("P.Value", "logFC", "pct.1", "pct.2", 
                                   "adj.P.Val", "Tag", "geneSymbol")
            my.markers$Test <- paste0(my.markers$Tag,"-Rest")
            my.markers$AveExpr <- NA
            my.markers$t <- NA
            my.markers$B <- NA
            my.markers$Assay <- my_assay
            my.markers$Method <- 'Wilcox'
            my.markers <- my.markers[,c("logFC", "AveExpr", "t", 
                                        "P.Value", "adj.P.Val", "B", "Test","Tag",
                                        "geneSymbol",'Method','Assay',"pct.1","pct.2")]
            si_cell_markers_tmp <- rbind(si_cell_markers_tmp,my.markers)
          }
        }
        rm(my.markers)
      }
      rm(mycolm)
    }else if(si_compute_cellmarker ==F & isTruthy(si_cell_markers)){
      mycols <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Test","Tag","geneSymbol")
      cf <- NULL
      ext <- tools::file_ext(si_cell_markers)
      if(ext=='txt'){
        cf = read.delim(si_cell_markers,as.is = T,colClasses = "#",sep="\t")
      }else if(ext == 'csv'){
        cf = read.csv2(si_cell_markers,as.is = T,colClasses = "#")
      }
      if(sum(mycols%in%names(cf))!=length(mycols)){
        warning("One or more columns are missing in the provided cell marker file. Please ensure that correct file is provided. Skipping this file and moving on!")
      }else{
        mycols <- append(mycols,
                         names(cf)[!names(cf)%in%mycols]
        )
        si_cell_markers_tmp <- cf[,mycols]
      }
      rm(cf,ext)
    }
    if(isTRUE(nrow(si_cell_markers_tmp)>0)){
      writedb(tableName=paste0(StudyName,"_Marker"),datafrm=si_cell_markers_tmp,db_address=db_address,OUTDIR=OUTDIR,overwrite=T)
    }
    rm(si_cell_markers_tmp)
    gc()
    plan("sequential")

  }
  gc()

  if(intractv==T){
    scProg$set(message = "Done!", value = 1)
  }
  message("Done!")
  
  return(1)
}


write.scstudy2.sqlitedb <- function(so,
                                    db_address=NULL,
                                    StudyName=NULL,
                                    Celltype = NULL,
                                    Reduction_map = NULL,
                                    Donors_VariableName=NULL,
                                    Disease_VariableName=NULL,
                                    Marker_Calc = F,
                                    Marker_Covariates = NULL,
                                    Marker_Precomputed=NULL,
                                    DE_Calc=F,
                                    DE_Precomputed = NULL,
                                    StudyDescr=NULL,
                                    Organism=NULL,
                                    Disease=NULL,
                                    ShortDescr=NULL,
                                    Tissue=NULL,
                                    PMID=NULL,
                                    GEO=NULL,
                                    StudyStatus='Internal',
                                    StudyRating='High',
                                    Continuous_Vars=NULL,
                                    Categorical_Vars=NULL,
                                    OUTDIR=NULL,
                                    signature_file=NULL,
                                    intractv=F){
  
  ##---- dependencies
  require(dplyr) || stop("Can't find dplyr package. Please install and restart! \n")
  require(Matrix) || stop("Can't find Matrix package. Please install and restart! \n")
  require(data.table) || stop("Can't find data.table package. Please install and restart! \n")
  require(RSQLite) || stop("Can't find RSQLite package. Please install and restart! \n")
  require(Seurat) || stop("Can't find Seurat package. Please install and restart! \n")
  
  if(intractv==T){
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = "assessing dependancies..", value = 0)
  }
  
  ##-- prechecks
  if(is.null(db_address)){
    stop('Please provide valid database address where the data is to be written')
  }
  if(length(grep("[.]db",db_address))==0){
    #db_address <- gsub("_sc$","",db_address)
    #db_address <- paste0(db_address,"_sc.db")
    db_address <- paste0(db_address,".db")
  }else{
    #db_address <- gsub(".db","_sc.db",db_address)
  }
  message("database address: ",db_address,"\n")
  if(is.null(StudyName)){
    stop('Please provide StudyName using which you would like to access this particualr data in the SingleCellApp')
  }
  if(is.null(Organism)){
    stop('Please provide Organism name for recording the study.')
  }
  if(is.null(Disease)){
    stop('Please provide the disease term for recording the study.')
  }
  if(is.null(ShortDescr)){
    stop('Please provide short one-line description for recording the study.')
  }
  if(is.null(StudyDescr)){
    warning('No study description is provided. You will see this section as blank in the Application\n')  
  }else{
    message(StudyDescr,"\n\n")
  }
  if(is.null(Donors_VariableName)){
    warning('Donors variable name is not provided. This information in the summary section will remain blank\n')  
  }else{
    message('Identifier for donors: ',Donors_VariableName,"\n\n")
  }
  if(is.null(Tissue)){
    warning('Tissue used for the experiment is not provided. This information in the summary section will remain blank\n')  
  }else{
    message('Identifier for tissue: ',Tissue,"\n\n")
  }
  if(is.null(Disease_VariableName)){
    warning('Disease variable (e.g. experimental condition that is not control/healthy etc.) for the experiment is not provided. This information in the summary section will remain blank\n')  
  }else{
    message('Identifier for disease: ',Disease_VariableName,"\n\n")
  }
  if(is.null(PMID)){
    message('(If public), the PMID for the experiment is not provided. This information in the summary section will remain blank\n')  
  }else{
    message('PMID: ',PMID,"\n\n")
  }
  if(is.null(GEO)){
    message('(If public), the data link for the experiment is not provided. This information in the summary section will remain blank\n')  
  }else{
    message('Data link (if public): ',GEO,"\n\n")
  }
  if(is.null(StudyStatus)){
    warning('If public/Internal, the status for the experiment is not provided. This information in the summary section will remain blank\n')  
  }else{
    message('Study status: ',StudyStatus,"\n\n")
  }
  if(is.null(StudyRating)){
    warning('If public/Internal, the status rating (e.g. high confidence/ poor study etc.) for the experiment is not provided. This information in the summary section will remain blank\n')  
  }else{
    message('Study rating: ',StudyRating,"\n\n")
  }
  if(is.null(Continuous_Vars)){
    warning('No continuous varibles in the metadata file are identified. This information in the summary section will remain blank\n')  
  }else{
    qq <- unlist(strsplit(as.character(Continuous_Vars),","))
    if(sum(qq%in%names(so@meta.data))!=length(qq) ){
      stop(" could not find provided continuous variables in the metadata file\n") 
    }
    message('Continuous variables in the metadata: ',Continuous_Vars,"\n\n")
  }
  if(is.null(Celltype)){
    stop("cell type column is not declared in the input. Please declare the column and rerun the function\n")
  }else{
    message(' cell type variable in the metadata: ', Celltype, "\n\n")
    so@meta.data <- local({
      meta <- so@meta.data
      meta$SVcell_type = meta[,Celltype]
      meta
    })
  }
  if(is.null(Categorical_Vars)){
    stop('No categorical varibles in the metadata file are identified. This information in the summary section will remain blank\n')  
  }else{
    qq <- unlist(strsplit(as.character(Categorical_Vars),","))
    if(sum(qq%in%names(so@meta.data))!=length(qq) ){
      stop(" could not find provided categorical variables in the metadata file\n") 
    }
    Categorical_Vars = paste0(qq,collapse = ",")
    message('Categorical variables in the metadata: ',Categorical_Vars,"\n\n")
  }
  if(!is.null(Marker_Covariates)){
    message(" Covariates for which to calculate marker genes: ", Marker_Covariates,"\n")
  }
  message("... 03.0: categorical variables before running functions: ",Categorical_Vars,"\n")
  
  study <- data.frame(Database=ifelse(!is.null(StudyName),as.character(StudyName),NA),
                      Description=ifelse(!is.null(StudyDescr),as.character(StudyDescr),NA),
                      SampleSize=ifelse(!is.na(Donors_VariableName),length(unique(so@meta.data[,Donors_VariableName])),NA),
                      TISSUES=ifelse(!is.null(Tissue),as.character(Tissue),NA),
                      DISEASE=ifelse(!is.null(Disease),as.character(Disease),NA),
                      ORGANISM=ifelse(!is.null(Organism),as.character(Organism),NA),
                      SHORTDESCR=ifelse(!is.null(ShortDescr),as.character(ShortDescr),NA),
                      DISEASE_Variable=ifelse(!is.null(Disease_VariableName),as.character(Disease_VariableName),NA), 
                      DONOR_Variable=ifelse(!is.na(Donors_VariableName),as.character(Donors_VariableName),NA), 
                      PMID=ifelse(!is.null(PMID),as.character(PMID),NA),
                      GEO=ifelse(!is.null(GEO),as.character(GEO),NA),
                      STATUS=ifelse(!is.null(StudyStatus),as.character(StudyStatus),NA),
                      RATING=ifelse(!is.null(StudyRating),as.character(StudyRating),NA),
                      CONTVAR=ifelse(!is.null(Continuous_Vars),as.character(Continuous_Vars),NA),
                      CATVAR=ifelse(!is.null(Categorical_Vars),as.character(Categorical_Vars),NA))
  
  message("... 03.1: categorical variables before running functions: ",study$CATVAR,"\n")
  
  
  gc()
  if(intractv==T){
    scProg$set(message = "performing analyses..", value = 0.2)
  }
  
  seurat2sqlite(so=so,
                si_study=study,
                si_reduction=Reduction_map,
                si_compute_cellmarker = Marker_Calc,
                si_cell_markers=Marker_Precomputed,
                si_cell_col = Marker_Covariates,
                si_compute_diseasemarker = DE_Calc,
                si_disease_markers=DE_Precomputed,
                si_disease_col=Disease_VariableName,
                si_celltype=Celltype,#'SVcell_type',
                si_donor=Donors_VariableName,
                db_address=db_address,
                OUTDIR=OUTDIR,
                signature_file=signature_file,
                intractv = intractv)
  
  if(intractv==T){
    scProg$set(message = 'analyses completed and SQLite file written!', value = 1)
  }
  return(1)
}



readInpscRNAseq <- function(filepath=NULL, intractv=F){
  
  ext <- tools::file_ext(filepath)
  #VARS$outdir <- dirname(filepath)
  message("observed file extention: ",ext,"\n")
  validate(need(ext %in% c("h5ad","rds"), "Please upload a correct file!"))
  #output$direxists <- renderPrint({ cat('outdir: ',VARS$outdir) })
  so = NULL
  if(ext=='rds'){
    if(intractv==T){
      scProg <- shiny::Progress$new()
      on.exit(scProg$close())
      scProg$set(message = "reading the file..", value = 0)
    }
    so <- local({
      #-- this part simplifies the Seurat object and reduces size, removes version incompatibility etc
      sa = readRDS(filepath)
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
      if(length(scvis)>1){
        message('Multiple data reductions found. Will populate them.\n')
        if(intractv==T){
          scProg$set(message = "'Multiple data reductions found. Will populate them.", value = 0.7)
        }
      } 
      for(i in 1:length(scvis)){
        message(" ..adding coordinates: ",gsub("^X_","",names(scvis)[i]),"\n")
        so[[gsub("^X_","",names(scvis)[i])]] = CreateDimReducObject(embeddings = scvis[[i]]@cell.embeddings, key = gsub("^X_","",names(scvis)[i]))
      }
      so
    })
    #variables$inptfile <- paste0(' --seurat \"',unname(inFile$datapath),"\"")
    if(intractv==T){
      scProg$set(message = "done!", value = 1)
    }
  }else if(ext =='h5ad'){
    if(intractv==T){
      scProg <- shiny::Progress$new()
      on.exit(scProg$close())
      scProg$set(message = "reading the file..", value = 0)
    }
    so <- local({
      ad <- anndata::read_h5ad(filename = filepath)
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
      if(length(scvis)>1){
        message('Multiple data reductions found. Will populate them.\n')
        if(intractv==T){
          scProg$set(message = "'Multiple data reductions found. Will populate them.", value = 0.7)
        }
      } 
      for(i in 1:length(scvis)){
        message(" ..adding coordinates: ",gsub("^X_","",names(scvis)[i]),"\n")
        so[[gsub("^X_","",names(scvis)[i])]] = CreateDimReducObject(embeddings = scvis[[i]], key = gsub("^X_","",names(scvis)[i]))
      }
      so
    })
    #variables$inptfile <- paste0(' --h5ad \"',unname(inFile$datapath),"\"")
    if(intractv==T){
      scProg$set(message = "done!", value = 1)
    }
  }
  return(so)
}




