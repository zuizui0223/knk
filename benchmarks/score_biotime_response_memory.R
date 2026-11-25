#!/usr/bin/env Rscript
library(data.table)
args <- commandArgs(TRUE)
stopifnot(length(args)==3L)
dir.create(args[3], recursive=TRUE, showWarnings=FALSE)
norm <- function(x) { s <- sum(x); if(s>0) x/s else rep(0,length(x)) }
tv <- function(a,b) .5*sum(abs(a-b))
top <- function(x) if(max(x)>0) which(x==max(x)) else integer()
score <- function(z) {
  # z: four years x two sample halves x species; all counts are mean per event.
  A <- z[,1,]; B <- z[,2,]
  d <- B[4,]-B[3,]; dA <- A[4,]-A[3,]; g <- sum(abs(d))
  h <- abs(A[2,]-A[1,]); last <- A[3,]; avg <- colMeans(A[1:3,,drop=FALSE])
  base <- list(norm(last),norm(avg),norm(sqrt(avg)))
  p <- norm(h)
  pred <- c(base,list(if(sum(h)>0) (base[[2]]+p)/2 else base[[2]],
                      if(sum(h)>0) (base[[3]]+p)/2 else base[[3]]))
  # With empty predictors, zero vector has TV 0.5 but is not a distribution: treat as undefined.
  loss <- vapply(pred,function(v) if(g>0 && sum(v)>0) tv(norm(abs(d)),v) else NA_real_,0.)
  tt <- top(last); tr <- top(abs(d))
  grossA <- sum(abs(dA))
  ca <- if(g>0) 1-abs(sum(d))/g else NA_real_
  sg <- (sum(sign(dA)*d)+sum(sign(d)*dA))/2
  sn <- (sign(sum(dA))*sum(d)+sign(sum(d))*sum(dA))/2
  pd <- norm(B[3,]); qd <- norm(B[4,])
  out <- c(M0=loss[1],M1=loss[2],M2=loss[3],M3=loss[4],M4=loss[5],
    current_top_response_share=if(g>0 && length(tt)>0) mean(abs(d)[tt])/g else NA_real_,
    current_top_overlap=if(length(tr)>0 && length(tt)>0) length(intersect(tt,tr))/(length(tt)*length(tr)) else NA_real_,
    lag_top_overlap=if(length(tr)>0 && length(top(h))>0) length(intersect(top(h),tr))/(length(top(h))*length(tr)) else NA_real_,
    gross=g,grossA=grossA,cancellation=ca,signed_gross=sg,signed_net=sn,signed_gap=sg-sn,
    namedTV=if(sum(pd)>0&&sum(qd)>0) tv(pd,qd) else NA_real_,
    rankTV=if(sum(pd)>0&&sum(qd)>0) tv(sort(pd),sort(qd)) else NA_real_,
    lag_defined=as.numeric(sum(h)>0), species=length(d))
  stopifnot(all(is.na(loss)|loss>=0 & loss<=1+1e-9))
  as.list(out)
}
# Synthetic: swap identity, perfect cancellation, no information leakage.
z <- array(0,c(4,2,2))
for(j in 1:2) z[,j,] <- rbind(c(5,0),c(5,5),c(10,0),c(5,5))
s <- score(z); stopifnot(abs(s$cancellation-1)<1e-10, s$gross==10, s$rankTV<=s$namedTV+1e-12)
z[4,1,] <- c(500,500); s2 <- score(z)
stopifnot(identical(s[c("M0","M1","M2","M3","M4")],s2[c("M0","M1","M2","M3","M4")]))
cat("Synthetic metric checks passed\n")

meta <- fread(args[2],encoding="UTF-8")
names(meta) <- toupper(names(meta))
needm <- c("STUDY_ID","ABUNDANCE_TYPE","REALM")
if(!all(needm %in% names(meta))) stop(paste("metadata schema:",paste(names(meta),collapse=" ")))
meta[, STUDY_ID:=as.character(STUDY_ID)]
meta[, type:=tolower(trimws(as.character(ABUNDANCE_TYPE)))]
# Exclude eBird-labelled records wherever explicit; do not assume other birds are eBird.
charcols <- names(meta)[vapply(meta,is.character,TRUE)]
flag <- rep(FALSE,nrow(meta))
for (col in charcols) flag <- flag | grepl("\\bebird\\b",meta[[col]],ignore.case=TRUE)
meta[, eb:=flag]
eligible <- meta[,.(count_ok=all(!is.na(type)&type=="count"),eb=any(eb),realm=paste(sort(unique(na.omit(REALM))),collapse="|")),by=STUDY_ID][count_ok & !eb & !grepl("\\|",realm)]
fwrite(eligible,file.path(args[3],"metadata_eligible.csv"))
raw <- readRDS(args[1]); setDT(raw)
names(raw) <- toupper(names(raw))
need <- c("ID_ALL_RAW_DATA","STUDY_ID","SAMPLE_DESC","LATITUDE","LONGITUDE","YEAR","MONTH","ABUNDANCE","VALID_NAME","RESOLUTION")
if(!all(need %in% names(raw))) stop(paste("raw schema:",paste(names(raw),collapse=" ")))
if(anyDuplicated(raw$ID_ALL_RAW_DATA)) stop("duplicate raw record IDs")
if(!"DEPTH" %in% names(raw)) raw[,DEPTH:=NA_real_]
if(!"PLOT" %in% names(raw)) raw[,PLOT:=NA_character_]
audit <- data.table(stage="raw",rows=nrow(raw),studies=uniqueN(raw$STUDY_ID))
raw[,STUDY_ID:=as.character(STUDY_ID)]
raw <- raw[STUDY_ID %chin% eligible$STUDY_ID]
audit <- rbind(audit,data.table(stage="count_metadata",rows=nrow(raw),studies=uniqueN(raw$STUDY_ID)))
raw <- raw[,.(STUDY_ID,SAMPLE_DESC=as.character(SAMPLE_DESC),LATITUDE=as.numeric(LATITUDE),LONGITUDE=as.numeric(LONGITUDE),YEAR=as.integer(YEAR),MONTH=as.integer(MONTH),ABUNDANCE=as.numeric(ABUNDANCE),SPECIES=as.character(VALID_NAME),RESOLUTION=tolower(as.character(RESOLUTION)),DEPTH=as.character(DEPTH),PLOT=as.character(PLOT))]
raw <- raw[!is.na(SAMPLE_DESC) & nzchar(SAMPLE_DESC) & is.finite(LATITUDE)&abs(LATITUDE)<=90 & is.finite(LONGITUDE)&abs(LONGITUDE)<=180 & !is.na(YEAR)&MONTH %in% 1:12]
# Audit events with invalid count or non-species rows, do not silently drop mixed records.
raw[, good:=is.finite(ABUNDANCE)&ABUNDANCE>=0&!is.na(SPECIES)&nzchar(SPECIES)&!is.na(RESOLUTION)&RESOLUTION=="species"]
ev <- raw[,.(year=uniqueN(YEAR),month=uniqueN(MONTH),lat=uniqueN(LATITUDE),lon=uniqueN(LONGITUDE),depth=uniqueN(DEPTH),plot=uniqueN(PLOT),good=all(good)),by=.(STUDY_ID,SAMPLE_DESC)]
ev <- ev[year==1&month==1&lat==1&lon==1&depth==1&plot==1&good==TRUE,.(STUDY_ID,SAMPLE_DESC)]
raw <- raw[ev,on=.(STUDY_ID,SAMPLE_DESC),nomatch=0]
if(nrow(raw)==0) stop("no eligible events after strict audit")
audit <- rbind(audit,data.table(stage="strict_events",rows=nrow(raw),studies=uniqueN(raw$STUDY_ID)))
fwrite(audit,file.path(args[3],"source_audit.csv"))
rm(ev); gc()
# Encode study-stratum and source event; aggregate duplicated species within a sampling event.
raw[, exact:=paste(sprintf("%.8f",LATITUDE),sprintf("%.8f",LONGITUDE),DEPTH,PLOT,MONTH,sep="|")]
# Spherical equal-area 100km bins (not exact ellipsoidal EPSG6933).
raw[, grid:=paste(floor(6371*cos(pi/6)*LONGITUDE*pi/180/100),floor(6371*sin(LATITUDE*pi/180)/cos(pi/6)/100),MONTH,sep="|")]
set.seed(20261125)
run_mode <- function(mode) {
  events <- unique(raw[,.(STUDY_ID,SAMPLE_DESC,YEAR,site=get(mode))])
  # select exactly 4 events per year prior to any outcome scoring
  eligible_year <- events[, .N,by=.(STUDY_ID,site,YEAR)][N>=4]
  wins <- eligible_year[, {
    ys <- sort(YEAR); starts <- integer(); i <- 1L
    while(i<=length(ys)) {
      if(all((ys[i]+0:3) %in% ys)) {starts<-c(starts,ys[i]); i<-i+4L} else i<-i+1L
    }
    .(start=starts)
  },by=.(STUDY_ID,site)]
  if(nrow(wins)==0) return(data.table())
  years <- wins[,.(YEAR=start+0:3),by=.(STUDY_ID,site,start)]
  e <- events[years,on=.(STUDY_ID,site,YEAR),nomatch=0]
  setorder(e,STUDY_ID,site,YEAR,SAMPLE_DESC)
  sel <- e[, { ix<-sample.int(.N,4L); .(SAMPLE_DESC=SAMPLE_DESC[ix],half=rep(1:2,each=2L)) },by=.(STUDY_ID,site,start,YEAR)]
  stopifnot(nrow(sel)==nrow(wins)*16L)
  # Only selected events, using the audited unique study/event key.
  agg <- raw[sel,on=.(STUDY_ID,SAMPLE_DESC),nomatch=0][,.(value=sum(ABUNDANCE)/2),by=.(STUDY_ID,site,start,i=YEAR-start+1L,half,SPECIES)]
  # Each distinct species is aligned over four years and two independent halves.
  results <- agg[, {
    taxa<-sort(unique(SPECIES)); arr<-array(0,c(4,2,length(taxa)))
    arr[cbind(i,half,match(SPECIES,taxa))]<-value
    # drop=FALSE in score not needed once at least two taxa present.
    if(length(taxa)<2) NULL else score(arr)
  },by=.(STUDY_ID,site,start)]
  if(nrow(results)>0) results[, mode:=mode]
  results
}
ans <- rbindlist(list(run_mode("exact"),run_mode("grid")),fill=TRUE)
if(nrow(ans)==0) stop("no eligible four-year windows with >=2 species")
ans <- eligible[ans,on="STUDY_ID",nomatch=0]
ans[,gain_primary:=M2-M4]; ans[,gain_sensitivity:=M1-M3]
fwrite(ans,file.path(args[3],"window_metrics.csv"))
# Paired losses require valid predictions and nonzero target; count all undefined windows.
study <- ans[,.(windows=.N,scored=sum(is.finite(gain_primary)),lag_defined=mean(lag_defined),
   M2=mean(M2,na.rm=TRUE),M4=mean(M4,na.rm=TRUE),gain=mean(gain_primary,na.rm=TRUE),
   gain_sens=mean(gain_sensitivity,na.rm=TRUE),cancellation=mean(cancellation,na.rm=TRUE),
   signed_gap=mean(signed_gap,na.rm=TRUE),top_overlap=mean(current_top_overlap,na.rm=TRUE),
   lag_top_overlap=mean(lag_top_overlap,na.rm=TRUE),
   driver_share=mean(current_top_response_share,na.rm=TRUE)),by=.(mode,realm,STUDY_ID)]
fwrite(study,file.path(args[3],"study_metrics.csv"))
set.seed(5511)
summary <- study[, {
  g<-gain[is.finite(gain)]
  ci<-if(length(g)>=10) quantile(replicate(2000,mean(sample(g,length(g),replace=TRUE))),c(.025,.975)) else c(NA,NA)
  .(n_studies=.N,scored_studies=length(g),n_windows=sum(windows),scored_windows=sum(scored),
    gain_mean=mean(g),study_improved=sum(g>0),CI_low=ci[1],CI_high=ci[2],
    loo_min=if(length(g)>1) min((sum(g)-g)/(length(g)-1)) else NA_real_,
    loo_max=if(length(g)>1) max((sum(g)-g)/(length(g)-1)) else NA_real_,
    mean_cancellation=mean(cancellation,na.rm=TRUE),mean_signed_gap=mean(signed_gap,na.rm=TRUE),
    mean_top_overlap=mean(top_overlap,na.rm=TRUE),mean_lag_top_overlap=mean(lag_top_overlap,na.rm=TRUE))
},by=mode]
fwrite(summary,file.path(args[3],"summary.csv"))
print(audit);print(summary);print(study[,.(n_studies=.N,n_windows=sum(windows),n_scored=sum(scored)),by=.(mode,realm)])
