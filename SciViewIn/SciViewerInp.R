#!/usr/bin/env Rscript
require("optparse")
require("getopt")
source("utilities.R",local = TRUE,encoding = "UTF-8")
source("vars.R",local = TRUE,encoding = "UTF-8")

option_list = list(
  ##------- Input data files
  make_option(c("--h5ad"), type="character", default=NULL, 
              help="If the data file is in the H5AD format, then the file path [default= NULL]", metavar="character"),
  make_option(c("--seurat"), type="character", default=NULL, 
              help="If the data file is in the Seurat format, then the file path [default= NULL]", metavar="character"),
  
  ##------- Mandatory inputs
  make_option(c("--db_address"), type="character", default=NULL, 
              help="(*) full address where the database file is to be written [default= NULL]", metavar="character"),
  make_option(c("--StudyName"), type="character", default=NULL, 
              help="(*) name using which you would like to identify the study in the application [default= NULL]", metavar="character"),
  make_option(c("--celltypeCol"), type="character", default=NULL, 
              help="(*) name of cell type annotation column in the metadata file [default= NULL]", metavar="character"),
  make_option(c("--org"), type="character", default=NULL, 
              help="(*) organism name [default= NULL]", metavar="character"),
  make_option(c("--diseaseName"), type="character", default=NULL, 
              help="(*) disease term [default= NULL]", metavar="character"),
  make_option(c("--shotdescr"), type="character", default=NULL, 
              help="(*) short one-liner description of the study [default= NULL]", metavar="character"),
  
  ##------- Optional inputs from section-2
  make_option(c("--reduction"), type="character", default=NULL, 
              help="(*) Default coordinate system to display.  [default= NULL]", metavar="character"),
  make_option(c("--donorCol"), type="character", default=NULL, 
              help="Name of donor annotation column in the metadata file [default= NULL]", metavar="character"),
  
  ##------- Optional inputs from section-3
  make_option(c("--Continuous_Vars"), type="character", default=NULL, 
              help="continuous variables identified from the metadata file, this info will be used for data display [default= %default]", metavar="character"),
  make_option(c("--Categorical_Vars"), type="character", default=NULL, 
              help="categorical variables identified from the metadata file, this info will be used for data display. column name 'cell_type' is required [default= %default]", metavar="character"),
  
  ##------- Inputs for DE/Markers calculation, section-4
  make_option(c("--compute_markers"), type="logical", default=T, 
              help="whether to compute markers for cell type and any additional provided covariate [default= %default]", metavar="logical"),
  make_option(c("--covariates"), type="character", default=NULL, 
              help="Comma-separated list of covariates from the metadata file for which markers genes to be calculated [default= NULL]", metavar="character"),
  make_option(c("--marker_file"), type="character", default=NULL, 
              help="Precomputed and formatted text file path containing marker information [default= NULL]", metavar="character"),
  make_option(c("--compute_deg"), type="logical", default=F, 
              help="whether to compute differential gene expression [default= %default]", metavar="logical"),
  make_option(c("--disease_var"), type="character", default=NULL, 
              help="Column header for the test variable [default= NULL]", metavar="character"),
  make_option(c("--deg_file"), type="character", default=NULL, 
              help="Precomputed and formatted text file path containing differentially expressed genes between conditions [default= NULL]", metavar="character"),
  
  ##------- Study details, section-5.1
  make_option(c("--StudyDescr"), type="character", default=NULL, 
              help="detailed description of the study, this info will be reflected in the summary section [default= %default]", metavar="character"),
  make_option(c("--Tissue"), type="character", default=NULL, 
              help="tissue used in the study, this info will be reflected in the summary section", metavar="character"),
  make_option(c("--PMID"), type="character", default=NULL, 
              help="publication id if available, this info will be reflected in the summary section [default= %default]", metavar="character"),
  make_option(c("--GEO"), type="character", default=NULL, 
              help="data repository link if available, this info will be reflected in the summary section [default= %default]", metavar="character"),
  make_option(c("--StudyStatus"), type="character", default='Internal', 
              help="status of the study (e.g. public/internal), this info will be reflected in the summary section [default= %default]", metavar="character"),
  make_option(c("--StudyRating"), type="character", default='High', 
              help="your rating for the study (e.g. high/poor/average), this info will be reflected in the summary section [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


###--------------------- H5AD to data objects
FILEPATH = NULL
if(!is.null(opt$h5ad)){
  FILEPATH=opt$h5ad
}else if(!is.null(opt$seurat)){
  FILEPATH=opt$seurat
}
OUTDIR = dirname(FILEPATH)

###--------------------- Input tests
if(is.null(FILEPATH)){
  stop(' provided input format is not supported!')
}
if(is.null(opt$db_address)){
  stop('path to write database is not provided!')
}
if(is.null(opt$celltypeCol)){
  stop('column name for cell type annotations not provided!')
}
if(is.null(opt$StudyName)){
  stop('name of the study to be used for writing the tables in the database is not specified!')
}

###--------------------- run
if(!is.null(FILEPATH)){
  so <- readInpscRNAseq(filepath=FILEPATH,intractv=F)
}

# precomputed marker file path
Marker_Precomputed = NULL
if(!is.null(opt$marker_file)){
  Marker_Precomputed = opt$marker_file
}

# precomputed de file path
DE_Precomputed = NULL
if(!is.null(opt$deg_file)){
  DE_Precomputed = opt$deg_file
}

##--- signature file
signature_file=NULL
if(file.exists("data/msigdb_signatures.txt")){
  signature_file <- read.delim("data/msigdb_signatures.txt",header=F,stringsAsFactors = F)
  if(sum(c("V1", "V2", "V3")%in%names(signature_file))==3){
    message(" ..gene signature file exits, will be used for computing the signatures\n")
  }else{
    message(" ..gene signature file exits, but not correctly formatted. Skipping signature calculations.\n")
  }
}else{
  message(" ..gene signature file does not exits. Skipping signature calculations.\n")
}



###--------------------- generate database
write.scstudy2.sqlitedb(so = so,
                        db_address=opt$db_address,
                        StudyName=opt$StudyName,
                        Donors_VariableName = opt$donorCol,
                        Celltype = opt$celltypeCol,
                        Reduction_map = opt$reduction,
                        Organism=opt$org,
                        Disease=opt$diseaseName,
                        ShortDescr=opt$shotdescr,
                        
                        Marker_Calc = opt$compute_markers,
                        Marker_Covariates = opt$covariates,
                        Marker_Precomputed=Marker_Precomputed,
                        
                        DE_Calc = opt$compute_deg,
                        Disease_VariableName=opt$disease_var,
                        DE_Precomputed=DE_Precomputed,
                        
                        StudyDescr=opt$StudyDescr,
                        Tissue=opt$Tissue,
                        PMID=opt$PMID,
                        GEO=opt$GEO,
                        StudyStatus=opt$StudyStatus,
                        StudyRating=opt$StudyRating,
                        
                        Continuous_Vars=opt$Continuous_Vars,
                        Categorical_Vars=opt$Categorical_Vars,
                        
                        OUTDIR = OUTDIR,
                        signature_file=signature_file)






