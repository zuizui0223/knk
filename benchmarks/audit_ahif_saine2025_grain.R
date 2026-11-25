# Source data inspection only: no source author analysis code is executed.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
root <- args[[1]]
out <- args[[2]]
datafile <- file.path(root,"data_meta_otu.Rdata")
stopifnot(file.exists(datafile), unname(tools::md5sum(datafile))=="4dbc27f81ffdcf1925887c7e0805b031")
e <- new.env(parent=emptyenv())
load(datafile,envir=e)
cat("RDATA_OBJECTS\n")
for(n in ls(e)){
  x <- get(n,envir=e)
  cat("OBJECT",n,"CLASS",paste(class(x),collapse=","),"DIM",paste(dim(x),collapse="x"),"LENGTH",length(x),"\n")
  if(is.data.frame(x)) {
    cat("FIELDS",n,paste(names(x),collapse=" | "),"\n")
    if(n %in% c("meta","read.counts","taxonomy")){
      for(k in names(x)){
        z <- x[[k]]
        if(is.atomic(z) && length(unique(z))<=30) {
          cat("LEVELS",n,k,paste(names(table(z,useNA="ifany")),as.integer(table(z,useNA="ifany")),sep="=",collapse=" | "),"\n")
        }
      }
      print(utils::head(x,3))
    }
  }
}
stopifnot(all(c("meta","read.counts","otu.table.plausible","taxonomy") %in% ls(e)))
meta <- e$meta
X <- e$otu.table.plausible
rc <- e$read.counts
tax <- e$taxonomy
cat("SOURCE_RAW_SHAPE",nrow(meta),nrow(X),ncol(X),nrow(rc),nrow(tax),"\n")
cat("ROWNAME_ALIGNMENT",identical(rownames(X),rownames(meta)),identical(rownames(X),rownames(rc)),identical(colnames(X),rownames(tax)),"\n")
cat("MATRIX_COUNTS",sum(X>0),sum(is.na(X)),min(X,na.rm=TRUE),max(X,na.rm=TRUE),"\n")
stopifnot(nrow(meta)==1380,nrow(X)==1380,ncol(X)==4662,length(unique(meta$RunningLogID))==460)
stopifnot(nrow(rc)==1380, "filt_nread" %in% names(rc))
sel <- rc$filt_nread>=10000
cat("READ_QC",sum(sel),sum(!sel),"\n")
cat("META_TREATMENT",paste(names(table(meta$InocSpecies)),as.integer(table(meta$InocSpecies)),sep="=",collapse=" | "),"\n")
cat("AUTHOR_S2_FIRST_110_LINES\n")
s2 <- file.path(root,"S2_define_HMSC_models_and_export_to_Hmsc-HPC.R")
stopifnot(unname(tools::md5sum(s2))=="a3dfe8ae32340e470a31cfa7bac7f9dd")
cat(paste(sprintf("%d %s",1:110,readLines(s2,warn=FALSE)[1:110]),collapse="\n"),"\n")
writeLines(c("status=SOURCE_RDATA_GRANULARITY_AUDIT_ONLY",paste0("raw_rows=",nrow(X)),paste0("raw_otus=",ncol(X)),paste0("qc_rows=",sum(sel)),paste0("physical_logs=",length(unique(meta$RunningLogID)))),out)
