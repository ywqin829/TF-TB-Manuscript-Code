#!/usr/bin/env Rscript
# Preselected TF contrasts with uncertainty, using the exact composition PCs of script 16.
source("scripts/00_config.R")
suppressPackageStartupMessages({library(limma);library(dplyr);library(readr);library(tibble);library(dorothea)})
inputs <- readRDS(out_path("results","analysis_inputs.rds"))
activity <- list(GSE40553=readRDS(out_path("results","gse40553_discovery_results.rds"))$tf_activity,
                 GSE31348=readRDS(out_path("results","gse31348_development_results.rds"))$tf_activity,
                 GSE89403=readRDS(out_path("results","gse89403_outcome_results.rds"))$tf_activity)
meta <- list(GSE40553=inputs$gse40553$pheno,GSE31348=inputs$gse31348$metadata,GSE89403=inputs$gse89403$metadata)
make_pcs <- function(comp) {
  v <- log(pmin(pmax(comp,1e-4),1-1e-4)/(1-pmin(pmax(comp,1e-4),1-1e-4)))
  v <- apply(v,2,function(z){z[!is.finite(z)] <- median(z[is.finite(z)]);z})
  v <- v[,apply(v,2,var)>0,drop=FALSE]
  z <- prcomp(scale(v),center=TRUE,scale.=FALSE)$x[,1:3,drop=FALSE]
  colnames(z) <- paste0("composition_PC",1:3); z
}
fit_one <- function(ds, adjusted) {
  # The consensus correlation and empirical-Bayes prior must use the original full TF panel.
  y <- activity[[ds]]
  m <- meta[[ds]][colnames(y),,drop=FALSE]
  m$SubjectID <- factor(m$SubjectID)
  if(ds=="GSE40553") m$Group <- factor(m$Group,levels=c("LTB","TB_0M","TB_0.5M","TB_2M","TB_6M","TB_12M"))
  if(ds=="GSE31348") m$TimePoint <- factor(m$TimePoint,levels=c("Week 0","Week 1","Week 2","Week 4","Week 26"))
  if(ds=="GSE89403") {
    m$Group <- factor(m$Group,levels=c("Successful_Treatment","Poor_Outcome"))
    m$TimePoint <- factor(m$TimePoint,levels=c("0_Baseline","1_Week1","2_Week4","3_Week24"))
    m$GroupTime <- interaction(m$Group,m$TimePoint,drop=TRUE)
  }
  term <- switch(ds,GSE40553="Group",GSE31348="TimePoint",GSE89403="GroupTime")
  if(adjusted) {
    c <- read.delim(out_path("results",paste0(tolower(ds),"_quantiseq.tsv")),check.names=FALSE)
    rownames(c) <- c$SampleID;c$SampleID <- NULL
    m <- cbind(m,make_pcs(as.matrix(c[rownames(m),,drop=FALSE])))
  }
  f <- as.formula(paste("~0+",term,if(adjusted) "+composition_PC1+composition_PC2+composition_PC3" else ""))
  d <- model.matrix(f,data=m);colnames(d) <- make.names(colnames(d))
  cm <- matrix(0,nrow=ncol(d),ncol=1,dimnames=list(colnames(d),"contrast"))
  if(ds=="GSE40553") cm[c("GroupTB_6M","GroupTB_0M"),1] <- c(1,-1)
  if(ds=="GSE31348") cm[c("TimePointWeek.26","TimePointWeek.0"),1] <- c(1,-1)
  if(ds=="GSE89403") {
    poor <- grep("^GroupTimePoor_Outcome\\.",rownames(cm),value=TRUE)
    good <- grep("^GroupTimeSuccessful_Treatment\\.",rownames(cm),value=TRUE)
    cm[poor,1] <- 1/length(poor);cm[good,1] <- -1/length(good)
  }
  dc <- duplicateCorrelation(y,d,block=m$SubjectID)
  fit <- eBayes(contrasts.fit(lmFit(y,d,block=m$SubjectID,correlation=dc$consensus.correlation),cm))
  tibble(Dataset=ds,TF=rownames(y),Model=if(adjusted) "composition_adjusted" else "unadjusted",
    Estimate=as.numeric(fit$coefficients[,1]),SE=as.numeric(fit$stdev.unscaled[,1]*fit$sigma),
    DF=as.numeric(fit$df.total),P_value=as.numeric(fit$p.value[,1])) %>%
    mutate(CI_low=Estimate-qt(.975,DF)*SE,CI_high=Estimate+qt(.975,DF)*SE) %>%
    filter(TF %in% c("PAX5","TBX21","MAF"))
}
z <- bind_rows(lapply(names(activity),function(ds) bind_rows(fit_one(ds,FALSE),fit_one(ds,TRUE))))
write_tsv(z,out_path("results","key_tf_composition_effects.tsv"))
# Assess direct target/marker overlap for PAX5; the marker score in script 34 excludes CD79A.
reg <- dorothea::dorothea_hs
targets <- unique(reg$target[reg$tf=="PAX5" & reg$confidence %in% c("A","B","C")])
markers <- c("CD79A","MS4A1","CD19","CD79B")
write_tsv(tibble(Marker=markers,In_PAX5_regulon=markers %in% targets,
                 In_independent_marker_score=markers %in% c("MS4A1","CD19","CD79B")),
          out_path("results","pax5_bcell_marker_target_overlap.tsv"))
print(z)
# Outcome contrast sensitivity to B-cell markers outside the PAX5 A/B/C regulon.
ds <- "GSE89403"
y <- activity[[ds]]
m <- meta[[ds]][colnames(y),,drop=FALSE]
m$SubjectID <- factor(m$SubjectID)
m$Group <- factor(m$Group,levels=c("Successful_Treatment","Poor_Outcome"))
m$TimePoint <- factor(m$TimePoint,levels=c("0_Baseline","1_Week1","2_Week4","3_Week24"))
m$GroupTime <- interaction(m$Group,m$TimePoint,drop=TRUE)
expr <- inputs$gse89403$expression[,colnames(y),drop=FALSE]
fit_marker <- function(marker_genes,label) {
  present <- intersect(marker_genes,rownames(expr))
  if(length(present)!=length(marker_genes)) stop("Missing marker genes")
  standardized <- t(scale(t(expr[present,,drop=FALSE])))
  m$marker_score <- as.numeric(colMeans(standardized))
  d <- model.matrix(~0+GroupTime+marker_score,data=m)
  colnames(d) <- make.names(colnames(d))
  cm <- matrix(0,nrow=ncol(d),ncol=1,dimnames=list(colnames(d),"contrast"))
  poor <- grep("^GroupTimePoor_Outcome\\.",rownames(cm),value=TRUE)
  good <- grep("^GroupTimeSuccessful_Treatment\\.",rownames(cm),value=TRUE)
  cm[poor,1] <- 1/length(poor); cm[good,1] <- -1/length(good)
  dc <- duplicateCorrelation(y,d,block=m$SubjectID)
  fit <- eBayes(contrasts.fit(lmFit(y,d,block=m$SubjectID,correlation=dc$consensus.correlation),cm))
  i <- match("PAX5",rownames(y)); estimate <- fit$coefficients[i,1]
  se <- fit$stdev.unscaled[i,1]*fit$sigma[i]; df <- fit$df.total[i]
  tibble(Model=label,Marker_genes=paste(present,collapse=";"),Estimate=estimate,
         CI_low=estimate-qt(.975,df)*se,CI_high=estimate+qt(.975,df)*se,
         P_value=fit$p.value[i,1])
}
marker_results <- bind_rows(
  fit_marker(c("CD79A","MS4A1","CD19","CD79B"),"four_markers"),
  fit_marker(c("MS4A1","CD79B"),"nonoverlap_markers"))
write_tsv(marker_results,out_path("results","pax5_nonoverlap_marker_sensitivity.tsv"))
print(marker_results)
