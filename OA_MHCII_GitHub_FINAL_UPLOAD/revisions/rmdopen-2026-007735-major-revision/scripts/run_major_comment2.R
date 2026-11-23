#!/usr/bin/env Rscript
set.seed(20261113)
cmd <- commandArgs(FALSE)
script_file <- sub('^--file=', '', cmd[grep('^--file=', cmd)][1])
root <- normalizePath(file.path(dirname(script_file), '..', '..'), mustWork=TRUE)
out <- file.path(root, 'revision_R1_major2')
dir.create(file.path(out,'results'), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out,'figures'), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out,'reports'), recursive=TRUE, showWarnings=FALSE)
inp <- read.csv(file.path(out,'results/canonical_matched_input.csv'), stringsAsFactors=FALSE, check.names=FALSE)
inp$group <- factor(inp$group, levels=c('Control','OA'))
inp$module <- as.numeric(inp$module_index_score)
inp$A <- as.numeric(inp$Reference_A_Myeloid_APC)
inp$Bv2 <- as.numeric(inp$Reference_Bv2_Myeloid_APC)
datasets <- c('GSE1919','GSE12021','GSE55235','GSE55457','GSE82107','GSE206848','GSE89408','GSE143514')
refs <- c('Reference_A','Reference_Bv2')
dir.create(file.path(out,'logs'), recursive=TRUE, showWarnings=FALSE)
log <- file(file.path(out,'logs/run.log'),'wt'); sink(log, type='output'); sink(log, type='message')
on.exit({sink(type='message'); sink(type='output'); close(log)}, add=TRUE)

spearman_boot <- function(x,y,B=2000) {
  ok <- is.finite(x)&is.finite(y); x<-x[ok]; y<-y[ok]; n<-length(x)
  if(n<3 || sd(x)==0 || sd(y)==0) return(c(rho=NA, p=NA, ci_low=NA, ci_high=NA, n=n, valid='NO_VARIATION_OR_SMALL_N'))
  rho <- suppressWarnings(cor(x,y,method='spearman'))
  p <- suppressWarnings(cor.test(x,y,method='spearman',exact=FALSE)$p.value)
  boot <- replicate(B, {ii <- sample.int(n,n,replace=TRUE); xx<-x[ii]; yy<-y[ii]; if(sd(xx)==0||sd(yy)==0) NA_real_ else suppressWarnings(cor(xx,yy,method='spearman'))})
  boot<-boot[is.finite(boot)]
  ci<-if(length(boot)>=100) quantile(boot,c(.025,.975),names=FALSE,type=7) else c(NA,NA)
  c(rho=rho,p=p,ci_low=ci[1],ci_high=ci[2],n=n,valid=if(length(boot)>=100)'YES' else 'BOOTSTRAP_UNSTABLE')
}
partial_spearman <- function(x,y,g) {
  rx <- rank(x,ties.method='average'); ry <- rank(y,ties.method='average'); z <- as.numeric(g=='OA')
  if(sd(rx)==0 || sd(ry)==0 || sd(z)==0) return(list(rho=NA,p=NA,n=length(x),valid='NO_VARIATION'))
  ex <- resid(lm(rx~z)); ey <- resid(lm(ry~z)); rho <- cor(ex,ey)
  df <- length(x)-3
  p <- if(df>0 && abs(rho)<1) 2*pt(-abs(rho)*sqrt(df/(1-rho^2)),df=df) else NA_real_
  list(rho=rho,p=p,n=length(x),valid='YES')
}
partial_boot <- function(x,y,g,B=2000) {
  base <- partial_spearman(x,y,g); if(!is.finite(base$rho)) return(list(rho=NA,p=NA,n=length(x),valid=base$valid,ci_low=NA,ci_high=NA,bootstrap_valid=0))
  idx0<-which(g=='Control'); idx1<-which(g=='OA')
  boot<-replicate(B,{ii<-c(sample(idx0,length(idx0),replace=TRUE),sample(idx1,length(idx1),replace=TRUE)); z<-partial_spearman(x[ii],y[ii],g[ii]); z$rho})
  boot<-boot[is.finite(boot)]; ci<-if(length(boot)>=100) quantile(boot,c(.025,.975),names=FALSE,type=7) else c(NA,NA)
  list(rho=base$rho,p=base$p,n=base$n,valid=base$valid,ci_low=ci[1],ci_high=ci[2],bootstrap_valid=length(boot))
}
make_row <- function(ds, ref, analysis, sub, x,y) {
  z<-if(analysis=='group_specific') spearman_boot(x,y) else partial_boot(x,y,sub)
  data.frame(dataset=ds, reference=ref, analysis=analysis, subgroup=if(analysis=='group_specific')as.character(sub[1]) else 'OA_adjusted',
             n=as.integer(if(is.list(z))z$n else z['n']), rho=as.numeric(if(is.list(z))z$rho else z['rho']), p_value=as.numeric(if(is.list(z))z$p else z['p']), ci_low=as.numeric(if(is.list(z))z$ci_low else z['ci_low']), ci_high=as.numeric(if(is.list(z))z$ci_high else z['ci_high']),
             valid=as.character(if(is.list(z))z$valid else z['valid']), bootstrap_valid=ifelse(is.list(z),as.numeric(z$bootstrap_valid),NA), stringsAsFactors=FALSE)
}
rows<-list(); k<-1
for(ds in datasets) for(ref in refs) {
  d<-inp[inp$dataset==ds,]; y<-d[[if(ref=='Reference_A')'A' else 'Bv2']]
  for(gp in c('OA','Control')) { ss<-d$group==gp; rows[[k]]<-make_row(ds,ref,'group_specific',d$group[ss],d$module[ss],y[ss]); k<-k+1 }
  rows[[k]]<-make_row(ds,ref,'OA_adjusted',d$group,d$module,y); k<-k+1
}
cohort <- do.call(rbind,rows)
cohort$meta_eligible <- with(cohort, n>=6 & valid=='YES')
write.csv(cohort,file.path(out,'results/cohort_level_correlations.csv'),row.names=FALSE,quote=TRUE)

library(metafor)
meta_one <- function(subset, analysis, ref, variance_rule, subgroup_filter=NULL) {
  q<-subset[subset$analysis==analysis & subset$reference==ref & is.finite(subset$rho) & subset$meta_eligible,]
  if(!is.null(subgroup_filter)) q<-q[q$subgroup==subgroup_filter,]
  if(nrow(q)<3) return(data.frame(analysis=analysis,reference=ref,subgroup=ifelse(is.null(subgroup_filter),'OA_adjusted',subgroup_filter),k=nrow(q),pooled_rho=NA,ci_low=NA,ci_high=NA,pi_low=NA,pi_high=NA,p_value=NA,HK_SE=NA,tau2=NA,I2=NA,status='NOT_ESTIMATED'))
  yi<-atanh(pmax(pmin(q$rho,0.999999),-0.999999)); vi<-if(analysis=='OA_adjusted')1/(q$n-4) else 1/(q$n-3)
  fit<-rma(yi=yi,vi=vi,method='REML',test='knha')
  pr<-predict(fit,transf=tanh)
  data.frame(analysis=analysis,reference=ref,subgroup=ifelse(is.null(subgroup_filter),'OA_adjusted',subgroup_filter),k=nrow(q),pooled_rho=tanh(as.numeric(coef(fit))),ci_low=tanh(fit$ci.lb),ci_high=tanh(fit$ci.ub),pi_low=as.numeric(pr$pi.lb),pi_high=as.numeric(pr$pi.ub),p_value=fit$pval,HK_SE=fit$se,tau2=fit$tau2,I2=fit$I2,status=paste0('REML_HK_',variance_rule),stringsAsFactors=FALSE)
}
metas<-do.call(rbind,c(lapply(refs,function(r)meta_one(cohort,'group_specific',r,'1/(n-3)','OA')),lapply(refs,function(r)meta_one(cohort,'group_specific',r,'1/(n-3)','Control')),lapply(refs,function(r)meta_one(cohort,'OA_adjusted',r,'1/(n-4)'))))
write.csv(metas,file.path(out,'results/group_specific_meta_results.csv'),row.names=FALSE,quote=TRUE)

loo<-list(); z<-1
for(an in c('group_specific','OA_adjusted')) for(ref in refs) for(sg in if(an=='group_specific') c('OA','Control') else 'OA_adjusted') {
  q<-cohort[cohort$analysis==an & cohort$reference==ref & cohort$subgroup==sg & cohort$meta_eligible,]
  for(ds in q$dataset) { qq<-q[q$dataset!=ds,]; if(nrow(qq)>=3) { yi<-atanh(qq$rho); vi<-if(an=='OA_adjusted')1/(qq$n-4) else 1/(qq$n-3); fit<-rma(yi=yi,vi=vi,method='REML',test='knha'); pr<-predict(fit,transf=tanh); loo[[z]]<-data.frame(analysis=an,reference=ref,subgroup=sg,omitted_dataset=ds,k=nrow(qq),pooled_rho=tanh(coef(fit)),ci_low=tanh(fit$ci.lb),ci_high=tanh(fit$ci.ub),pi_low=pr$pi.lb,pi_high=pr$pi.ub); z<-z+1} }
}
write.csv(do.call(rbind,loo),file.path(out,'results/leave_one_dataset_out.csv'),row.names=FALSE,quote=TRUE)

orig_rows<-list(); oi<-1
for(ds in datasets) { dd<-inp[inp$dataset==ds,]; for(ref in refs) { yy<-dd[[if(ref=='Reference_A')'A' else 'Bv2']]; rr<-suppressWarnings(cor(dd$module,yy,method='spearman')); orig_rows[[oi]]<-data.frame(dataset=ds,reference=ref,analysis='original_combined',subgroup='OA+Control',n=nrow(dd),rho=rr,p_value=NA,ci_low=NA,ci_high=NA,valid='YES',bootstrap_valid=NA,meta_eligible=TRUE); oi<-oi+1 } }
orig_rows<-do.call(rbind,orig_rows)
orig <- data.frame(reference=c('Reference_A','Reference_Bv2'), expected=c(0.572213432117994,0.548270616956164), stringsAsFactors=FALSE)
orig$reproduced <- sapply(orig$reference,function(r) meta_one(orig_rows,'original_combined',r,'1/(n-3)')$pooled_rho)
orig$absolute_difference <- abs(orig$expected-orig$reproduced)
orig$within_tolerance <- orig$absolute_difference<1e-8
write.csv(orig,file.path(out,'results/original_pooled_reproduction.csv'),row.names=FALSE,quote=TRUE)

summary_rows <- rbind(
  data.frame(check='records',value=as.character(nrow(inp))),
  data.frame(check='unique_dataset_sample_keys',value=as.character(length(unique(paste(inp$dataset,inp$sample_id))))),
  data.frame(check='finite_module',value=as.character(all(is.finite(inp$module)))),
  data.frame(check='finite_A',value=as.character(all(is.finite(inp$A)))),
  data.frame(check='finite_Bv2',value=as.character(all(is.finite(inp$Bv2)))),
  data.frame(check='fraction_A_bounded',value=as.character(all(inp$A>=0&inp$A<=1))),
  data.frame(check='fraction_Bv2_bounded',value=as.character(all(inp$Bv2>=0&inp$Bv2<=1))),
  data.frame(check='OA_control_label_consistent',value=as.character(all(inp$group==ifelse(inp$OA==1,'OA','Control')))),
  data.frame(check='original_pooled_Reference_A_reproduced',value=as.character(orig$within_tolerance[orig$reference=='Reference_A'])),
  data.frame(check='original_pooled_Reference_Bv2_reproduced',value=as.character(orig$within_tolerance[orig$reference=='Reference_Bv2'])),
  data.frame(check='group_specific_meta_k',value=paste(unique(metas$k[metas$analysis=='group_specific']),collapse=';')),
  data.frame(check='OA_adjusted_meta_k',value=paste(unique(metas$k[metas$analysis=='OA_adjusted']),collapse=';')),
  data.frame(check='bootstrap_replicates',value='2000; seed=20261113'),
  data.frame(check='partial_variance_approximation',value='1/(n-4)')
)
write.csv(summary_rows,file.path(out,'results/statistical_sensitivity_summary.csv'),row.names=FALSE,quote=TRUE)
saveRDS(list(cohort=cohort,meta=metas,session=sessionInfo()),file.path(out,'logs/analysis_objects.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'logs/sessionInfo.txt'))
sink(type='message'); sink(type='output'); close(log)
