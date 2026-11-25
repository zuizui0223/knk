#!/usr/bin/env Rscript
library(data.table)
a<-commandArgs(TRUE); stopifnot(length(a)==2L)
dir.create(a[2],recursive=TRUE,showWarnings=FALSE)
e<-new.env(parent=baseenv())
stopifnot(identical(load(file.path(a[1],"BTsubset_data.RData"),envir=e),"BTsubset_data"))
stopifnot(identical(load(file.path(a[1],"BTsubset_meta.RData"),envir=e),"BTsubset_meta"))
d<-as.data.table(e$BTsubset_data); m<-as.data.table(e$BTsubset_meta)
setnames(d,toupper(names(d))); setnames(m,toupper(names(m)))
stopifnot(nrow(d)==81084L,nrow(m)==12L,!anyDuplicated(d$ID_ALL_RAW_DATA))
d[,STUDY_ID:=as.character(STUDY_ID)]; m[,STUDY_ID:=as.character(STUDY_ID)]
stopifnot(all(unique(d$STUDY_ID)%chin%m$STUDY_ID))
# Conservative source-level bird exclusion; teaching metadata do not carry full citations.
avian<-function(x) grepl("bird|avian|\\baves\\b|ebird",as.character(x),ignore.case=TRUE)
excluded<-unique(c(m[avian(TAXA)|avian(TITLE),STUDY_ID],d[avian(TAXON),STUDY_ID]))
audit<-d[,.(rows=.N,taxa=paste(sort(unique(TAXON)),collapse="|")),by=STUDY_ID]
audit<-merge(audit,m[,.(STUDY_ID,REALM,ABUNDANCE_TYPE,TITLE)],by="STUDY_ID",all.x=TRUE)
audit[,excluded_avian:=STUDY_ID%chin%excluded]
fwrite(audit,file.path(a[2],"subset_source_studies.csv"))
print(audit)
saveRDS(d[!STUDY_ID%chin%excluded],file.path(a[1],"subset.rds"))
fwrite(m[!STUDY_ID%chin%excluded],file.path(a[1],"subset.csv"))
# Two nonoverlapping four-year windows in one site. Fixed 4 events/year; no count noise.
fix<-CJ(YEAR=2001:2008,rep=1:4,SPECIES=c("Species A","Species B"))
fix[,phase:=(YEAR-2001L)%%4L+1L]
fix[,ABUNDANCE:=ifelse(SPECIES=="Species A",c(5,5,10,5)[phase],c(0,5,0,5)[phase])]
fix[,SAMPLE_DESC:=paste(YEAR,rep,sep="_")]
fix[,ID_ALL_RAW_DATA:=.I]
fix[,":="(STUDY_ID="fixture",LATITUDE=0,LONGITUDE=0,MONTH=6L,VALID_NAME=SPECIES,RESOLUTION="species",DEPTH=0,PLOT="test")]
saveRDS(fix,file.path(a[1],"fixture.rds"))
fwrite(data.table(STUDY_ID="fixture",ABUNDANCE_TYPE="Count",REALM="Terrestrial"),file.path(a[1],"fixture.csv"))
