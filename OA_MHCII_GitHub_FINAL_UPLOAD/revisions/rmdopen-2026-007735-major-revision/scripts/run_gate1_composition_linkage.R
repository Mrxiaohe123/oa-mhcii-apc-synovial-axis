#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20261004)
suppressPackageStartupMessages({library(ggplot2); library(metafor)})

root <- "."
base <- file.path(root, "全部分析内容")
out <- file.path(base, "gate1_composition_linkage")
dir.create(file.path(out, "figures"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out, "source_data"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out, "logs"), recursive=TRUE, showWarnings=FALSE)

primary <- c("GSE1919","GSE12021","GSE55235","GSE55457","GSE82107","GSE206848","GSE89408","GSE143514")
all_ds <- c(primary, "GSE283079")
module <- c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","CD74","HLA-DMA","HLA-DMB","CIITA")
proxy <- c("LST1","TYROBP","FCER1G","AIF1","C1QA","C1QB","C1QC","CSF1R","CD68","FCGR3A","ITGAM","LYZ","SPI1","MS4A7","CTSS")
expr_dir <- file.path(base, "phase4f_reference_based_deconvolution_sensitivity/03_deconvolution_inputs")
frac_dir <- file.path(base, "phase4f_reference_based_deconvolution_sensitivity/05_music")

norm_names <- function(x) toupper(gsub("\\s+", "", trimws(x)))
read_expr <- function(ds) {
  f <- file.path(expr_dir, paste0(ds, "_bulk_expression_for_music.tsv"))
  x <- read.delim(f, check.names=FALSE, stringsAsFactors=FALSE)
  names(x)[1] <- "gene_symbol"
  x$gene_symbol <- norm_names(x$gene_symbol)
  x <- x[nzchar(x$gene_symbol),,drop=FALSE]
  num <- as.data.frame(lapply(x[-1], function(z) as.numeric(as.character(z))), check.names=FALSE)
  rownames(num) <- x$gene_symbol
  # Duplicate probes are averaged; this is a deterministic within-dataset operation.
  if (anyDuplicated(rownames(num))) {
    num <- aggregate(num, by=list(gene_symbol=rownames(num)), FUN=function(z) mean(z, na.rm=TRUE))
    rownames(num) <- num$gene_symbol; num$gene_symbol <- NULL
  }
  num
}
zscore <- function(x) { x <- as.numeric(x); if (sd(x, na.rm=TRUE)==0) rep(0,length(x)) else as.numeric(scale(x)) }
score_genes <- function(expr, genes) {
  g <- intersect(genes, rownames(expr)); if (!length(g)) return(list(score=rep(NA_real_, ncol(expr)), genes=character()))
  # sapply returns samples x genes because each gene vector is column-wise; average across genes per sample.
  list(score=rowMeans(sapply(g, function(z) zscore(expr[z,])), na.rm=TRUE), genes=g)
}
read_meta <- function(ds) {
  f <- file.path(expr_dir, paste0(ds, "_metadata_for_deconvolution.csv"))
  m <- read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
  m$sample <- as.character(m$sample)
  m$disease_group <- tolower(as.character(m$disease_group))
  m
}
read_frac <- function(ds) {
  f <- file.path(frac_dir, paste0(ds, "_MuSiC_cell_fractions.csv"))
  read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
}
ci_lm <- function(fit, term) {
  cc <- summary(fit)$coefficients
  if (!(term %in% rownames(cc))) return(c(NA,NA,NA,NA))
  b <- cc[term,"Estimate"]; se <- cc[term,"Std. Error"]; df <- fit$df.residual
  c(beta=b, lo=b-qt(.975,df)*se, hi=b+qt(.975,df)*se, p=cc[term,"Pr(>|t|)"])
}
boot_rho <- function(x,y,B=2000,seed=1) {
  set.seed(seed); n <- length(x); vals <- rep(NA_real_,B)
  for (i in seq_len(B)) { ii <- sample.int(n,n,replace=TRUE); if (length(unique(x[ii]))>1 && length(unique(y[ii]))>1) vals[i] <- suppressWarnings(cor(x[ii],y[ii],method="spearman")) }
  quantile(vals, c(.025,.975), na.rm=TRUE, names=FALSE)
}
safe_meta <- function(yi, sei, label, primary=TRUE) {
  ok <- is.finite(yi)&is.finite(sei)&sei>0; yi<-yi[ok]; sei<-sei[ok]
  if (length(yi)<3) return(data.frame(analysis=label,k=length(yi),pooled_effect=NA,ci_low=NA,ci_high=NA,p_value=NA,tau2=NA,I2=NA,Q=NA,pred_low=NA,pred_high=NA))
  fit <- rma(yi=yi, sei=sei, method="REML", test="knha")
  pr <- predict(fit, transf=if (grepl("rho",label)) tanh else identity)
  transf <- if (grepl("rho",label)) tanh(fit$b) else as.numeric(fit$b)
  ci <- if (grepl("rho",label)) tanh(c(fit$ci.lb,fit$ci.ub)) else c(fit$ci.lb,fit$ci.ub)
  data.frame(analysis=label,k=length(yi),pooled_effect=transf,ci_low=ci[1],ci_high=ci[2],p_value=fit$pval,tau2=fit$tau2,I2=fit$I2,Q=fit$QE,pred_low=pr$pi.lb,pred_high=pr$pi.ub)
}

rows <- list(); score_rows <- list(); att_rows <- list(); marker_rows <- list(); multi_rows <- list()
for (j in seq_along(all_ds)) {
  ds <- all_ds[j]; ex <- read_expr(ds); me <- read_meta(ds); fr <- read_frac(ds)
  ms <- score_genes(ex,module); ps <- score_genes(ex,proxy)
  names(ms$score) <- colnames(ex); names(ps$score) <- colnames(ex)
  dat <- data.frame(sample=colnames(ex), module=as.numeric(ms$score), proxy=as.numeric(ps$score))
  dat <- merge(dat, me[,c("sample","disease_group")], by="sample")
  fr$sample <- as.character(fr$sample); dat <- merge(dat,fr,by="sample")
  dat$OA <- as.integer(dat$disease_group=="oa")
  dat$dataset <- ds
  dat$role <- ifelse(ds %in% primary,"primary","external_holdout")
  dat <- dat[is.finite(dat$module)&is.finite(dat$Myeloid_APC),]
  write.csv(dat, file.path(out,"source_data",paste0(ds,"_sample_level_gate1_data.csv")),row.names=FALSE,na="")
  if (nrow(dat)>=4 && length(unique(dat$Myeloid_APC))>1 && length(unique(dat$module))>1) {
    ct <- suppressWarnings(cor.test(dat$module,dat$Myeloid_APC,method="spearman",exact=FALSE))
    bc <- boot_rho(dat$module,dat$Myeloid_APC,2000,20261004+j)
    fit <- lm(scale(module) ~ scale(Myeloid_APC), data=dat); lc <- ci_lm(fit,"scale(Myeloid_APC)")
    rows[[length(rows)+1]] <- data.frame(dataset=ds,role=dat$role[1],n=nrow(dat),n_OA=sum(dat$OA),n_control=sum(!dat$OA),module_genes_available=length(ms$genes),spearman_rho=unname(ct$estimate),spearman_ci_low=bc[1],spearman_ci_high=bc[2],spearman_p=ct$p.value,linear_beta=lc[1],linear_ci_low=lc[2],linear_ci_high=lc[3],linear_p=lc[4],r_squared=summary(fit)$r.squared,myeloid_min=min(dat$Myeloid_APC),myeloid_max=max(dat$Myeloid_APC),pathological_proportion=any(dat$Myeloid_APC<0|dat$Myeloid_APC>1))
    # Disease-state models on the locked module score.
    m1 <- lm(scale(module)~OA,data=dat); m2 <- lm(scale(module)~scale(Myeloid_APC),data=dat); m3 <- lm(scale(module)~OA+scale(Myeloid_APC),data=dat)
    a1<-ci_lm(m1,"OA"); a3<-ci_lm(m3,"OA")
    att_rows[[length(att_rows)+1]] <- data.frame(dataset=ds,role=dat$role[1],n=nrow(dat),model1_beta=a1[1],model1_low=a1[2],model1_high=a1[3],model1_p=a1[4],model2_myeloid_beta=ci_lm(m2,"scale(Myeloid_APC)")[1],model3_beta=a3[1],model3_low=a3[2],model3_high=a3[3],model3_p=a3[4],delta=a3[1]-a1[1],ratio=ifelse(abs(a1[1])>.1,a3[1]/a1[1],NA),attenuation_percent=ifelse(abs(a1[1])>.1,100*(1-a3[1]/a1[1]),NA),classification=ifelse(abs(a3[1])<abs(a1[1]),"attenuated",ifelse(abs(a3[1])>abs(a1[1]),"amplified","unchanged")))
    # Marker-proxy sensitivity uses the frozen non-overlapping proxy, separately from MuSiC abundance.
    if (length(ps$genes)>=2 && length(unique(dat$proxy))>1 && all(is.finite(dat$proxy))) {
      q1<-lm(scale(module)~OA,data=dat); q3<-lm(scale(module)~OA+scale(proxy),data=dat); b1<-ci_lm(q1,"OA"); b3<-ci_lm(q3,"OA")
      marker_rows[[length(marker_rows)+1]] <- data.frame(dataset=ds,role=dat$role[1],n=nrow(dat),proxy_genes_available=length(ps$genes),model1_beta=b1[1],model1_low=b1[2],model1_high=b1[3],model3_beta=b3[1],model3_low=b3[2],model3_high=b3[3],delta=b3[1]-b1[1],ratio=ifelse(abs(b1[1])>.1,b3[1]/b1[1],NA),attenuation_percent=ifelse(abs(b1[1])>.1,100*(1-b3[1]/b1[1]),NA),classification=ifelse(abs(b3[1])<abs(b1[1]),"attenuated",ifelse(abs(b3[1])>abs(b1[1]),"amplified","unchanged")))
    }
    comps <- intersect(c("Fibroblast","Endothelial","T_NK","B_cell"),names(dat)); if (length(comps)>=2) {
      dat_mc <- dat[complete.cases(dat[,c("module","OA","Myeloid_APC",comps),drop=FALSE]),,drop=FALSE]
      valid_mc <- nrow(dat_mc)>=5 && length(unique(dat_mc$OA))>=2 && sd(dat_mc$module)>0 && sd(dat_mc$Myeloid_APC)>0 && all(vapply(dat_mc[,c("Myeloid_APC",comps)], function(z) sd(z)>0, logical(1)))
      if (valid_mc) {
        fml <- as.formula(paste("scale(module)~OA+scale(Myeloid_APC)+",paste(sprintf("scale(%s)",comps),collapse="+")))
        mf <- tryCatch(lm(fml,data=dat_mc), error=function(e) NULL)
        if (!is.null(mf)) {
          vv <- tryCatch(max(diag(solve(cor(model.matrix(mf)[,-1,drop=FALSE])))),error=function(e)NA)
          multi_rows[[length(multi_rows)+1]] <- data.frame(dataset=ds,n=nrow(dat_mc),n_covariates=length(comps)+1,OA_beta=ci_lm(mf,"OA")[1],OA_low=ci_lm(mf,"OA")[2],OA_high=ci_lm(mf,"OA")[3],OA_p=ci_lm(mf,"OA")[4],max_VIF=vv,r_squared=summary(mf)$r.squared)
        }
      }
    }
    score_rows[[length(score_rows)+1]] <- dat
  }
}

direct <- do.call(rbind,rows); att <- do.call(rbind,att_rows); marker <- do.call(rbind,marker_rows); multi <- do.call(rbind,multi_rows); sample_all <- do.call(rbind,score_rows)
write.csv(direct,file.path(out,"source_data/module_vs_myeloid_dataset_results.csv"),row.names=FALSE)
write.csv(att,file.path(out,"source_data/oa_attenuation_dataset_results.csv"),row.names=FALSE)
write.csv(marker,file.path(out,"source_data/marker_proxy_sensitivity_dataset_results.csv"),row.names=FALSE)
write.csv(multi,file.path(out,"source_data/multicompartment_supportive_results.csv"),row.names=FALSE)

# OA-control composition contrast, calculated from MuSiC fractions with Hedges g.
comp <- do.call(rbind,lapply(split(sample_all,sample_all$dataset),function(d){
  oa<-d$Myeloid_APC[d$OA==1]; co<-d$Myeloid_APC[d$OA==0]; es<-escalc(measure="SMD",m1i=mean(oa),sd1i=sd(oa),n1i=length(oa),m2i=mean(co),sd2i=sd(co),n2i=length(co),append=TRUE)
  wt<-suppressWarnings(wilcox.test(oa,co,exact=FALSE)); data.frame(dataset=d$dataset[1],role=d$role[1],n_OA=length(oa),n_control=length(co),mean_OA=mean(oa),mean_control=mean(co),median_OA=median(oa),median_control=median(co),raw_difference=mean(oa)-mean(co),hedges_g=es$yi,hedges_se=sqrt(es$vi),ci_low=es$yi-qnorm(.975)*sqrt(es$vi),ci_high=es$yi+qnorm(.975)*sqrt(es$vi),wilcox_p=wt$p.value)
}))
write.csv(comp,file.path(out,"source_data/oa_vs_control_myeloid_results.csv"),row.names=FALSE)
elig <- data.frame(dataset=primary,
  OA_n=sapply(primary,function(z) sum(sample_all$dataset==z & sample_all$OA==1)),
  control_n=sapply(primary,function(z) sum(sample_all$dataset==z & sample_all$OA==0)),
  module_score_available=sapply(primary,function(z) z %in% direct$dataset),
  module_genes_available=sapply(primary,function(z) ifelse(z %in% direct$dataset,direct$module_genes_available[match(z,direct$dataset)],0)),
  marker_proxy_available=sapply(primary,function(z) z %in% marker$dataset),
  MuSiC_proportion_available=sapply(primary,function(z) z %in% comp$dataset),
  myeloid_APC_available=sapply(primary,function(z) z %in% comp$dataset),
  fibroblast_available=TRUE,endothelial_available=TRUE,T_NK_available=TRUE,B_cell_available=TRUE,
  exclusion_reason=ifelse(primary %in% direct$dataset,"","No complete module/Myeloid_APC sample-level data"))
write.csv(elig,file.path(out,"eligibility_table.csv"),row.names=FALSE)
write.csv(direct,file.path(out,"module_vs_myeloid_dataset_results.csv"),row.names=FALSE)
write.csv(comp,file.path(out,"oa_vs_control_myeloid_results.csv"),row.names=FALSE)
write.csv(att,file.path(out,"oa_attenuation_dataset_results.csv"),row.names=FALSE)
write.csv(marker,file.path(out,"marker_proxy_sensitivity_results.csv"),row.names=FALSE)

meta <- rbind(safe_meta(atanh(direct$spearman_rho[direct$role=="primary"]),1/sqrt(direct$n[direct$role=="primary"]-3),"primary_module_vs_myeloid_spearman"),safe_meta(comp$hedges_g[comp$role=="primary"],comp$hedges_se[comp$role=="primary"],"primary_OA_vs_control_myeloid_hedges_g"))
pool_lm <- function(x,label){x<-x[x$role=="primary"&is.finite(x$model1_beta),]; safe_meta(x$model1_beta,(x$model1_high-x$model1_low)/(2*qt(.975,x$n-2)),paste0(label,"_unadjusted"))}
pool_adj <- function(x,label){x<-x[x$role=="primary"&is.finite(x$model3_beta),]; safe_meta(x$model3_beta,(x$model3_high-x$model3_low)/(2*qt(.975,x$n-3)),paste0(label,"_adjusted"))}
meta <- rbind(meta,pool_lm(att,"primary_oa_module"),pool_adj(att,"primary_oa_module"),pool_lm(marker,"primary_marker_proxy"),pool_adj(marker,"primary_marker_proxy"))
write.csv(meta,file.path(out,"source_data/pooled_meta_results.csv"),row.names=FALSE)
write.csv(meta,file.path(out,"pooled_meta_results.csv"),row.names=FALSE)

# Leave-one-dataset-out direct-association sensitivity.
loo <- do.call(rbind,lapply(seq_len(nrow(direct)),function(i){d<-direct[-i,]; z<-safe_meta(atanh(d$spearman_rho[d$role=="primary"]),1/sqrt(d$n[d$role=="primary"]-3),paste0("LOO_",direct$dataset[i])); z$omitted_dataset<-direct$dataset[i];z}))
write.csv(loo,file.path(out,"source_data/leave_one_dataset_out_direct_meta.csv"),row.names=FALSE)

# Negative-control/QC inventory.
qc <- data.frame(check=c("MuSiC proportions within [0,1]","full module excluded from proxy","housekeeping comparator","primary dataset dominance check"),status=c(ifelse(all(sample_all$Myeloid_APC>=0&sample_all$Myeloid_APC<=1),"PASS","FAIL"),"PASS","NOT_AVAILABLE","SEE_LOO_TABLE"),detail=c("All observed Myeloid_APC fractions checked","Frozen 15-gene non-overlap proxy contains none of 10 module genes","No frozen housekeeping comparator was present in the audited inputs","Leave-one-dataset-out direct-association results exported"))
write.csv(qc,file.path(out,"source_data/negative_control_and_qc.csv"),row.names=FALSE)

# Figures G1-G5 (primary datasets only; holdout remains separately identifiable in source data).
sp <- sample_all[sample_all$role=="primary",]
p1 <- ggplot(sp,aes(Myeloid_APC,module,color=factor(OA,labels=c("Control","OA"))))+geom_point(size=1.5,alpha=.75)+geom_smooth(method="lm",se=FALSE,color="grey30",linewidth=.4)+facet_wrap(~dataset,scales="free")+scale_color_manual(values=c("grey55","#3B6E9E"))+theme_bw(base_size=9)+labs(x="MuSiC inferred Myeloid/APC proportion",y="10-gene module score",color=NULL,title="G1. Direct module–composition association")
ggsave(file.path(out,"figures/G1_module_vs_myeloid_scatter.png"),p1,width=9,height=6,dpi=300); ggsave(file.path(out,"figures/G1_module_vs_myeloid_scatter.pdf"),p1,width=9,height=6)
p2 <- ggplot(direct[direct$role=="primary",],aes(spearman_rho,reorder(dataset,spearman_rho)))+geom_vline(xintercept=0,lty=2,color="grey60") + geom_errorbarh(aes(xmin=spearman_ci_low,xmax=spearman_ci_high),height=.15,color="#3B6E9E")+geom_point(color="#3B6E9E",size=2)+theme_bw(base_size=9)+labs(x="Spearman rho (bootstrap 95% CI)",y=NULL,title="G2. Dataset-level direct association")
ggsave(file.path(out,"figures/G2_direct_association_forest.png"),p2,width=7,height=5,dpi=300); ggsave(file.path(out,"figures/G2_direct_association_forest.pdf"),p2,width=7,height=5)
p3 <- ggplot(comp[comp$role=="primary",],aes(hedges_g,reorder(dataset,hedges_g)))+geom_vline(xintercept=0,lty=2,color="grey60")+geom_errorbarh(aes(xmin=ci_low,xmax=ci_high),height=.15,color="#3B6E9E")+geom_point(color="#3B6E9E",size=2)+theme_bw(base_size=9)+labs(x="Hedges' g: OA minus control (95% CI)",y=NULL,title="G3. OA–control composition contrast")
ggsave(file.path(out,"figures/G3_oa_vs_control_composition.png"),p3,width=7,height=5,dpi=300); ggsave(file.path(out,"figures/G3_oa_vs_control_composition.pdf"),p3,width=7,height=5)
longa <- rbind(data.frame(dataset=att$dataset,estimate=att$model1_beta,low=att$model1_low,high=att$model1_high,model="Unadjusted"),data.frame(dataset=att$dataset,estimate=att$model3_beta,low=att$model3_low,high=att$model3_high,model="MuSiC-adjusted"))
p4 <- ggplot(longa,aes(model,estimate,group=dataset))+geom_hline(yintercept=0,lty=2,color="grey70")+geom_line(color="grey65",alpha=.7)+geom_point(aes(color=model),size=2)+geom_errorbar(aes(ymin=low,ymax=high),width=.08)+facet_wrap(~dataset,ncol=4)+scale_color_manual(values=c("grey45","#3B6E9E"))+theme_bw(base_size=9)+theme(legend.position="bottom")+labs(x=NULL,y="OA coefficient for standardized module (95% CI)",title="G4. OA coefficient before and after composition adjustment")
ggsave(file.path(out,"figures/G4_before_after_adjustment.png"),p4,width=9,height=6,dpi=300); ggsave(file.path(out,"figures/G4_before_after_adjustment.pdf"),p4,width=9,height=6)
sumdat <- data.frame(stage=c("OA–control composition","Module–composition association","Unadjusted OA→module coefficient","MuSiC-adjusted OA coefficient"),estimate=c(mean(comp$raw_difference[comp$role=="primary"]),mean(direct$spearman_rho[direct$role=="primary"]),mean(att$model1_beta[att$role=="primary"]),mean(att$model3_beta[att$role=="primary"])))
p5 <- ggplot(sumdat,aes(estimate,reorder(stage,estimate)))+geom_vline(xintercept=0,lty=2,color="grey70")+geom_point(size=3,color="#3B6E9E")+theme_bw(base_size=10)+labs(x="Descriptive summary estimate (unpooled mean)",y=NULL,title="G5. Gate 1 composition-sensitivity summary",subtitle="Observed associations; not a causal mediation diagram")
ggsave(file.path(out,"figures/G5_gate1_summary.png"),p5,width=8,height=4.5,dpi=300); ggsave(file.path(out,"figures/G5_gate1_summary.pdf"),p5,width=8,height=4.5)

writeLines(capture.output(sessionInfo()),file.path(out,"logs/sessionInfo.txt"))
writeLines(c(paste("Gate 1 run:",format(Sys.time())),"Seed: 20261004","Primary datasets: 8; external holdout: GSE283079","Existing manuscript and ART files were not modified."),file.path(out,"logs/run.log"))

cat("Gate 1 computation complete\n")
cat("direct datasets:",nrow(direct)," composition datasets:",nrow(comp)," attenuation datasets:",nrow(att),"\n")
cat("positive module-composition associations:",sum(direct$role=="primary" & direct$spearman_rho>0),"/",sum(direct$role=="primary"),"\n")
cat("higher OA composition means:",sum(comp$role=="primary" & comp$raw_difference>0),"/",sum(comp$role=="primary"),"\n")
cat("OA coefficient attenuation:",sum(att$role=="primary" & att$classification=="attenuated"),"/",sum(att$role=="primary"),"\n")
print(meta)
