# Audit and export only original verified objects. Temporary data stay in CI.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
root <- args[[1]]
out <- args[[2]]
datafile <- file.path(root,"data_meta_otu.Rdata")
s2 <- file.path(root,"S2_define_HMSC_models_and_export_to_Hmsc-HPC.R")
stopifnot(unname(tools::md5sum(datafile))=="4dbc27f81ffdcf1925887c7e0805b031",
          unname(tools::md5sum(s2))=="a3dfe8ae32340e470a31cfa7bac7f9dd")
cat("AUTHOR_S2_FIRST_110_LINES\n")
s2lines <- readLines(s2,warn=FALSE)
cat(paste(sprintf("%d %s",1:110,s2lines[1:110]),collapse="\n"),"\n")
e <- new.env(parent=emptyenv()); load(datafile,envir=e)
stopifnot(all(c("meta","read.counts.plausible","otu.table.plausible","taxonomy.plausible") %in% ls(e)))
meta <- e$meta; X <- as.matrix(e$otu.table.plausible); rc <- e$read.counts.plausible; tax <- e$taxonomy.plausible
cat("SOURCE_SHAPES",nrow(meta),dim(X),nrow(rc),nrow(tax),"\n")
stopifnot(identical(dim(X),c(1380L,4662L)),nrow(meta)==1380,nrow(rc)==1380,nrow(tax)==4662)
ids <- as.character(meta$SampleID); xr <- rownames(X); rr <- as.character(rc$sample); tr <- rownames(tax)
cat("KEY_HEADS",paste(head(ids,3),collapse="|"),paste(head(xr,3),collapse="|"),
    paste(head(rr,3),collapse="|"),paste(head(tr,3),collapse="|"),paste(head(colnames(X),3),collapse="|"),"\n")
stopifnot(!anyDuplicated(ids),!anyDuplicated(xr),!anyDuplicated(rr),!anyDuplicated(tr),
          !anyDuplicated(colnames(X)),setequal(ids,xr),setequal(ids,rr),setequal(colnames(X),tr))
X <- X[match(ids,xr),,drop=FALSE]; rc <- rc[match(ids,rr),,drop=FALSE]
tax <- tax[match(colnames(X),tr),,drop=FALSE]
stopifnot(identical(rownames(X),ids),identical(as.character(rc$sample),ids),
          identical(colnames(X),rownames(tax)),all(is.finite(X)),all(X>=0),all(X==floor(X)))
cat("READ_QC",sum(rc$filt_nread>=10000),sum(rc$filt_nread<10000),"\n")
stopifnot(all(is.finite(rc$filt_nread)),sum(rc$filt_nread>=10000)==1369)
# Author S2 lines 25--48, verified by the pinned source digest.
# These are deposited taxonomy labels, not guessed binomial matches.
inocs <- c("Antrodia_piceata","Antrodiella_citrinella","Fomitopsis_rosea",
           "Perenniporia_subacida","Physisporinus_crocatus","Postia_guttulata",
           "Skeletocutis_odora","Skeletocutis_stellae","Steccherinum_collabens")
spp <- c("Antrodia_piceata_813073","Flaviporus_citrinellus_106139",
         "Rhodofomes_roseus_127496","Perenniporia_subacida_335816",
         "Rigidoporus_crocatus_107233","Postia_guttulata_110917",
         "Skeletocutis_odora_106476","Skeletocutis_stellae_323595",
         "Junghuhnia_collabens_315979")
targets <- setdiff(unique(as.character(meta$InocSpecies)),"control")
stopifnot(setequal(targets,inocs))
matches <- lapply(spp,function(s) which(as.character(tax$species)==s))
cat("EXACT_AUTHOR_TARGET_MATCH_COUNTS",paste(inocs,lengths(matches),sep="=",collapse="|"),"\n")
stopifnot(all(lengths(matches)>0))
# S2's scalar assignment assumes one match per code; fail if ambiguity rather
# than silently picking first. This tests exactly the archived taxonomy.
stopifnot(all(lengths(matches)==1L))
ix <- unlist(matches,use.names=FALSE)
stopifnot(length(unique(ix))==9L)
mapped <- inocs[match(seq_len(ncol(X)),ix)]
sel <- seq_len(ncol(X)) %in% ix
cat("TARGET_OTU_MAP\n")
print(data.frame(otu=colnames(X)[ix],species=as.character(tax$species[ix]),mapped=inocs))
Y <- X[,!sel,drop=FALSE]
# Keep the fixed full roster before sample filtering; distinguish zero-support
# columns after read-depth QC from deliberately removed target OTUs.
qc <- rc$filt_nread>=10000
active <- colSums(Y[qc,,drop=FALSE])>0
cat("NONTARGET_ROSTER",ncol(Y),"QC_ACTIVE",sum(active),"QC_ZERO_COLUMNS",sum(!active),"\n")
cat("QC_ZERO_OTUS",paste(colnames(Y)[!active],collapse="|"),"\n")
m <- data.frame(sample=ids,year=as.character(meta$SamplingYear),site=as.character(meta$SiteNew),
                log=as.character(meta$RunningLogID),logtype=as.character(meta$LogType),
                treatment=as.character(meta$InocSpecies),reads=as.numeric(rc$filt_nread),
                stringsAsFactors=FALSE)
cat("LOGTYPE",paste(names(table(m$logtype)),as.integer(table(m$logtype)),sep="=",collapse="|"),"\n")
stopifnot(ncol(Y)==4653,length(unique(m$log))==460,
          setequal(unique(m$year),c("2019","2020","2021")),
          setequal(unique(m$logtype),c("B","R","U")))
cat("ZERO_NONTARGET_TOTAL",sum(rowSums(Y)==0),"\n")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(m,file.path(out,"meta.csv"),row.names=FALSE)
writeLines(colnames(Y),file.path(out,"otu_names.txt"))
writeLines(paste(colnames(X)[ix],inocs,spp,sep=":"),file.path(out,"excluded_target_otus.txt"))
con <- file(file.path(out,"counts.bin"),"wb")
writeBin(as.double(t(Y)),con,size=8,endian="little");close(con)
cat("GRAIN_EXPORT_OK",nrow(Y),ncol(Y),"\n")
