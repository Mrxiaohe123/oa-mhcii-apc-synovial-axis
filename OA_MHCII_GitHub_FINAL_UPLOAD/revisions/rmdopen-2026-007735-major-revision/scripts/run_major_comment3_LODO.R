#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE, warn=1, scipen=999)
suppressPackageStartupMessages({library(Matrix);library(SingleCellExperiment);library(MuSiC);library(metafor)})
cmd <- commandArgs(FALSE); script_file <- sub('^--file=','',cmd[grep('^--file=',cmd)][1])
ROOT <- normalizePath(file.path(dirname(script_file),'..','..'),mustWork=TRUE)
AN <- file.path(ROOT,'全部分析内容'); OUT <- file.path(ROOT,'revision_R1_major3')
dir.create(file.path(OUT,'results'),recursive=TRUE,showWarnings=FALSE); dir.create(file.path(OUT,'logs'),recursive=TRUE,showWarnings=FALSE)
LOG<-file(file.path(OUT,'logs/LODO_run.log'),'wt'); sink(LOG,split=TRUE); sink(LOG,type='message'); on.exit({sink(type='message');sink(type='output');close(LOG)},add=TRUE)
set.seed(20261113)
DATASETS<-c('GSE1919','GSE12021','GSE55235','GSE55457','GSE82107','GSE206848','GSE89408','GSE143514')
CLASSES<-c('Fibroblast','Myeloid_APC','B_cell','T_NK','Endothelial')
MODULE<-read.csv(file.path(ROOT,'revision_R1_major2/results/canonical_matched_input.csv'),check.names=FALSE,stringsAsFactors=FALSE)
MODULE$sample_id<-as.character(MODULE$sample_id); MODULE$OA<-as.integer(MODULE$OA); MODULE$module<-as.numeric(MODULE$module_index_score)

read_bulk <- function(ds) {
  if(ds!='GSE12021') { x<-read.table(file.path(AN,'reference_B_downstream_impact_v2/inputs/bulk',paste0(ds,'_bulk_expression_for_music.tsv')),header=TRUE,sep='\t',quote='',comment.char='',check.names=FALSE,row.names=1); x[is.na(x)]<-0; return(x) }
  # Use the frozen gene-level harmonized MuSiC input, then select the
  # corrected 19 GPL96 sample IDs. The raw GEO matrix is probe-level and
  # cannot be intersected directly with the gene-level reference.
  cs<-read.csv(file.path(AN,'gse12021_independent_unit_correction/GSE12021_GPL96_corrected_scores.csv'),check.names=FALSE,stringsAsFactors=FALSE)
  x<-read.table(file.path(AN,'phase4f_reference_based_deconvolution_sensitivity/03_deconvolution_inputs/GSE12021_bulk_expression_for_music.tsv'),header=TRUE,sep='\t',quote='',comment.char='',check.names=FALSE,row.names=1)
  keep<-as.character(cs$sample); stopifnot(all(keep%in%colnames(x)))
  x<-x[,keep,drop=FALSE]; x[is.na(x)]<-0; x
}

read_refA <- function(){readRDS(file.path(AN,'phase4f_reference_based_deconvolution_sensitivity/01_reference_construction/GSE248453_reference_primary_frozen10_excluded.rds'))}
read_refB <- function(){
  p<-file.path(AN,'reference_B_downstream_impact_v2/inputs/reference'); m<-readMM(file.path(p,'reference_B_counts.mtx')); g<-readLines(file.path(p,'reference_B_genes.txt')); md<-read.csv(file.path(p,'reference_B_cell_metadata.csv'),check.names=FALSE,stringsAsFactors=FALSE); rownames(m)<-g; colnames(m)<-md$cell_id; rownames(md)<-md$cell_id; SingleCellExperiment(assays=list(counts=m),colData=md)
}
fit_music <- function(expr,sce,ds,ref,omit_label,omit_donor) {
  tryCatch({
    if(!is.null(omit_donor)) sce<-sce[,as.character(colData(sce)$donor_id)!=omit_donor]
    common<-intersect(rownames(expr),rownames(sce)); expr<-expr[common,,drop=FALSE]; sce<-sce[common,]
    fit<-MuSiC::music_prop(bulk.mtx=as.matrix(expr),sc.sce=sce,clusters='coarse_class',samples='donor_id',verbose=FALSE)
    z<-as.matrix(fit$Est.prop.weighted); z<-z[,CLASSES,drop=FALSE]
    stopifnot(identical(rownames(z),colnames(expr)),all(is.finite(z)),all(z>=0&z<=1))
    r2<-as.numeric(fit$r.squared.full); stopifnot(length(r2)==nrow(z),all(is.finite(r2)))
    fr<-data.frame(reference=ref,omit_donor=omit_label,dataset=ds,sample_id=rownames(z),z,r_squared_full=r2,status='PASS',common_genes=length(common),stringsAsFactors=FALSE)
    diag<-data.frame(reference=ref,omit_donor=omit_label,dataset=ds,sample_id=rownames(z),fit_metric='MuSiC.r.squared.full',fit_value=r2,status='PASS',common_genes=length(common),stringsAsFactors=FALSE)
    list(fr=fr,diag=diag,fit=fit)
  },error=function(e){list(fr=NULL,diag=data.frame(reference=ref,omit_donor=omit_label,dataset=ds,sample_id=NA,fit_metric='MuSiC_r.squared.full',fit_value=NA,status=paste0('FAIL: ',conditionMessage(e)),common_genes=NA,stringsAsFactors=FALSE),fit=NULL)})
}

refs<-list(Reference_A=read_refA(),Reference_Bv2=read_refB())
donors<-lapply(refs,function(s)sort(unique(as.character(colData(s)$donor_id))))
all_fr<-list(); all_diag<-list(); failures<-list(); ii<-1
for(ref in names(refs)) for(omit in c('baseline',donors[[ref]])) {
  omit_donor<-if(omit=='baseline')NULL else omit
  cat('RUN',ref,omit,'\n')
  for(ds in DATASETS){
    z<-fit_music(read_bulk(ds),refs[[ref]],ds,ref,omit,omit_donor)
    if(!is.null(z$fr)){all_fr[[ii]]<-z$fr;ii<-ii+1}; all_diag[[length(all_diag)+1]]<-z$diag
    if(is.null(z$fr)) failures[[length(failures)+1]]<-z$diag
  }
}
fractions<-do.call(rbind,all_fr); diagnostics<-do.call(rbind,all_diag); write.csv(fractions,file.path(OUT,'results/LODO_sample_level_fractions.csv'),row.names=FALSE,quote=TRUE); write.csv(diagnostics,file.path(OUT,'results/deconvolution_fit_diagnostics.csv'),row.names=FALSE,quote=TRUE)
if(length(failures)) write.csv(do.call(rbind,failures),file.path(OUT,'results/LODO_failures.csv'),row.names=FALSE,quote=TRUE) else write.csv(data.frame(status='NO_FAILURES'),file.path(OUT,'results/LODO_failures.csv'),row.names=FALSE,quote=TRUE)

# Per-sample stability against the exact baseline run.
stab<-list(); si<-1
for(ref in names(refs)) for(omit in setdiff(unique(fractions$omit_donor),'baseline')) for(ds in DATASETS){
  b<-fractions[fractions$reference==ref&fractions$omit_donor=='baseline'&fractions$dataset==ds,]; q<-fractions[fractions$reference==ref&fractions$omit_donor==omit&fractions$dataset==ds,]; m<-merge(b,q,by=c('sample_id'),suffixes=c('_baseline','_lodo'))
  if(nrow(m)){d<-m$Myeloid_APC_lodo-m$Myeloid_APC_baseline; stab[[si]]<-data.frame(reference=ref,omit_donor=omit,dataset=ds,n=nrow(m),spearman_baseline_lodo=cor(m$Myeloid_APC_baseline,m$Myeloid_APC_lodo,method='spearman'),median_abs_difference=median(abs(d)),mean_abs_difference=mean(abs(d)),max_abs_difference=max(abs(d)),n_positive_delta=sum(d>0),stringsAsFactors=FALSE);si<-si+1}
}
write.csv(do.call(rbind,stab),file.path(OUT,'results/LODO_stability_summary.csv'),row.names=FALSE,quote=TRUE)

cohort_rho<-function(fr){
  out<-list();k<-1
  for(ref in names(refs)) for(omit in unique(fr$omit_donor[fr$reference==ref])) for(ds in DATASETS){
    x<-MODULE[MODULE$dataset==ds,]; y<-fr[fr$reference==ref&fr$omit_donor==omit&fr$dataset==ds,]; m<-merge(x,y,by.x='sample_id',by.y='sample_id'); r<-if(nrow(m)>=3) suppressWarnings(cor(m$module,m$Myeloid_APC,method='spearman')) else NA
    out[[k]]<-data.frame(reference=ref,omit_donor=omit,dataset=ds,n=nrow(m),spearman_rho=r,OA_n=sum(m$OA==1),Control_n=sum(m$OA==0),stringsAsFactors=FALSE);k<-k+1
  };do.call(rbind,out)
}
cr<-cohort_rho(fractions);write.csv(cr,file.path(OUT,'results/LODO_cohort_correlations.csv'),row.names=FALSE,quote=TRUE)
meta_rows<-list();mi<-1
for(ref in names(refs)) for(omit in unique(cr$omit_donor[cr$reference==ref])){
 q<-cr[cr$reference==ref&cr$omit_donor==omit&is.finite(cr$spearman_rho),]; yi<-atanh(pmax(pmin(q$spearman_rho,.999999),-.999999)); vi<-1/(q$n-3); fit<-rma(yi=yi,vi=vi,method='REML',test='knha'); pr<-predict(fit,transf=tanh); meta_rows[[mi]]<-data.frame(reference=ref,omit_donor=omit,k=nrow(q),pooled_rho=tanh(coef(fit)),ci_low=tanh(fit$ci.lb),ci_high=tanh(fit$ci.ub),pi_low=pr$pi.lb,pi_high=pr$pi.ub,p_value=fit$pval,tau2=fit$tau2,I2=fit$I2,status='REML_HK',stringsAsFactors=FALSE);mi<-mi+1
}
write.csv(do.call(rbind,meta_rows),file.path(OUT,'results/LODO_meta_results.csv'),row.names=FALSE,quote=TRUE)

reg_rows<-list();ri<-1
for(ref in names(refs)) for(omit in unique(fractions$omit_donor[fractions$reference==ref])) for(ds in DATASETS){
 x<-MODULE[MODULE$dataset==ds,]; y<-fractions[fractions$reference==ref&fractions$omit_donor==omit&fractions$dataset==ds,]; m<-merge(x,y,by='sample_id'); if(nrow(m)<4) next
 m$zmod<-as.numeric(scale(m$module)); m$zfrac<-as.numeric(scale(m$Myeloid_APC)); u<-lm(zmod~OA,data=m); a<-lm(zmod~OA+zfrac,data=m); bu<-coef(summary(u))['OA','Estimate']; seu<-coef(summary(u))['OA','Std. Error']; ba<-coef(summary(a))['OA','Estimate']; sea<-coef(summary(a))['OA','Std. Error']; reg_rows[[ri]]<-data.frame(reference=ref,omit_donor=omit,dataset=ds,n=nrow(m),unadjusted_beta=bu,unadjusted_se=seu,adjusted_beta=ba,adjusted_se=sea,attenuated=abs(ba)<abs(bu),direction_change=sign(ba)!=sign(bu),stringsAsFactors=FALSE);ri<-ri+1
}
reg<-do.call(rbind,reg_rows);write.csv(reg,file.path(OUT,'results/LODO_adjusted_regression.csv'),row.names=FALSE,quote=TRUE)
eff<-list();ei<-1
for(ref in names(refs)) for(omit in unique(reg$omit_donor[reg$reference==ref])){q<-reg[reg$reference==ref&reg$omit_donor==omit,]; ans<-list();for(an in c('unadjusted','adjusted')){yi<-if(an=='unadjusted')q$unadjusted_beta else q$adjusted_beta;vi<-(if(an=='unadjusted')q$unadjusted_se else q$adjusted_se)^2;fit<-rma(yi=yi,vi=vi,method='REML',test='knha');eff[[ei]]<-data.frame(reference=ref,omit_donor=omit,effect=an,k=nrow(q),pooled_beta=as.numeric(coef(fit)),ci_low=fit$ci.lb,ci_high=fit$ci.ub,p_value=fit$pval,tau2=fit$tau2,I2=fit$I2,attenuation_n=sum(q$attenuated),direction_change_n=sum(q$direction_change),status='REML_HK');ei<-ei+1}}
write.csv(do.call(rbind,eff),file.path(OUT,'results/LODO_adjusted_regression_meta.csv'),row.names=FALSE,quote=TRUE)
write.csv(data.frame(check=c('Reference_A_donors','Reference_Bv2_donors','LODO_expected_runs','LODO_successful_run_dataset_blocks','LODO_failures','baseline_fraction_rows'),value=c(length(donors$Reference_A),length(donors$Reference_Bv2),length(donors$Reference_A)+length(donors$Reference_Bv2),length(unique(paste(fractions$reference,fractions$omit_donor))),nrow(do.call(rbind,failures)),sum(fractions$omit_donor=='baseline'))),file.path(OUT,'results/analysis_audit.csv'),row.names=FALSE,quote=TRUE)
writeLines(capture.output(sessionInfo()),file.path(OUT,'logs/sessionInfo.txt'))
