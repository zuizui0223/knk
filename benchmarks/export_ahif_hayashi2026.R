# Read exact author rarefied outputs; do not reconstruct or redistribute sequencing reads.
args<-commandArgs(TRUE);stopifnot(length(args)==2L)
root<-args[1];out<-args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
meta<-read.csv(file.path(root,"expdata_rep.csv"),stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(meta)==4608L,!anyDuplicated(meta$Sample_name),setequal(unique(meta$Day),c(2,4,6,8,10,12)),setequal(meta$Source,c("soil","water")),setequal(meta$Inoculum_density,c(1,.1,.01,.001)))
meta$unit<-as.integer(sub("^.*_([0-9]+)$","\\1",meta$Sample_name))%%384L
stopifnot(!anyNA(meta$unit))
meta<-meta[meta$Day %in% c(2,4,6,8),]
stopifnot(nrow(meta)==3072L,all(table(meta$Source,meta$Inoculum_density,meta$Day)==96L),all(table(meta$Source,meta$unit)==4L))
for(g in split(meta,paste(meta$Source,meta$unit)))stopifnot(length(unique(g$Location))==1L,length(unique(g$Inoculum_density))==1L,length(unique(g$Plate_position))==1L)
stopifnot(all(meta$Plate_position==match(meta$Inoculum_density,c(1,.1,.01,.001))))
x<-readRDS(file.path(root,"List_seqtab_3000.rds"));y<-readRDS(file.path(root,"List_data_rep.rds"))
stopifnot(length(x)==7L,length(y)==8L)
audit<-list()
for(k in c(5,6)){
  tab<-as.matrix(x[[k]]); author<-y[[k]]
  stopifnot(!anyDuplicated(rownames(tab)),!anyDuplicated(colnames(tab)),all(is.finite(tab)),all(tab>=0),all(tab==floor(tab)),all(rowSums(tab)==3000))
  stopifnot(all(c("Sample_name","Day","Location","Source","Inoculum_density","Plate_position") %in% names(author)),!anyDuplicated(author$Sample_name))
  available<-meta[meta$Sample_name %in% rownames(tab),]
  u<-paste(available$Source,available$unit,sep=":")
  complete<-names(which(table(u)==4L));sub<-available[u %in% complete,]
  stopifnot(setequal(sub$Sample_name,author$Sample_name),nrow(sub)==4L*length(complete))
  author<-author[match(sub$Sample_name,author$Sample_name),,drop=FALSE]
  for(c in c("Day","Location","Source","Inoculum_density","Plate_position"))stopifnot(all(as.character(sub[[c]])==as.character(author[[c]])))
  obs<-tab[match(sub$Sample_name,rownames(tab)),,drop=FALSE]
  ac<-as.matrix(author[,7:ncol(author),drop=FALSE])
  stopifnot(identical(dim(obs),dim(ac)),all(obs==ac))
  fn<-if(k==6L)"otu99" else "asv"
  # Source table is exported only into ephemeral CI environment.
  write.csv(cbind(sub,as.data.frame(obs,check.names=FALSE)),file.path(out,paste0(fn,".csv")),row.names=FALSE,quote=TRUE)
  cat("HAYASHI_QC",fn,"all_table_rows",nrow(tab),"taxa",ncol(tab),"available_day2_8",nrow(available),"complete_cultures",length(complete),"complete_states",nrow(sub),"\n")
  print(with(sub,table(Source,Inoculum_density,Day)))
}
