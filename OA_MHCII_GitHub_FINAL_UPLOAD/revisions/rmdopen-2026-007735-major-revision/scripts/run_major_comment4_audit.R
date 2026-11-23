#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1,scipen=999)
suppressPackageStartupMessages(library(metafor))
cmd<-commandArgs(FALSE); sf<-sub('^--file=','',cmd[grep('^--file=',cmd)][1])
ROOT<-normalizePath(file.path(dirname(sf),'..','..'),mustWork=TRUE); OUT<-file.path(ROOT,'revision_R1_major4'); RES3<-file.path(ROOT,'revision_R1_major3/results'); A3<-file.path(ROOT,'revision_R1_major3/final_audit')
set.seed(20261113); dir.create(file.path(OUT,'results'),recursive=TRUE,showWarnings=FALSE); dir.create(file.path(OUT,'figures'),recursive=TRUE,showWarnings=FALSE)
DATASETS<-c('GSE1919','GSE12021','GSE55235','GSE55457','GSE82107','GSE206848','GSE89408','GSE143514')
reg<-read.csv(file.path(RES3,'LODO_adjusted_regression.csv'),stringsAsFactors=FALSE); base_reg<-subset(reg,omit_donor=='baseline')
meta<-read.csv(file.path(RES3,'LODO_adjusted_regression_meta.csv'),stringsAsFactors=FALSE); base_meta<-subset(meta,omit_donor=='baseline')
fractions<-subset(read.csv(file.path(RES3,'LODO_sample_level_fractions.csv'),stringsAsFactors=FALSE),omit_donor=='baseline')
mod<-read.csv(file.path(ROOT,'revision_R1_major2/results/canonical_matched_input.csv'),check.names=FALSE,stringsAsFactors=FALSE); mod$sample_id<-as.character(mod$sample_id); mod$OA<-as.integer(mod$OA); mod$module<-as.numeric(mod$module_index_score)

att<-base_reg[,c('reference','dataset','n','unadjusted_beta','unadjusted_se','adjusted_beta','adjusted_se','attenuated','direction_change')]
att$delta_abs<-abs(att$unadjusted_beta)-abs(att$adjusted_beta); att$delta_definition<-'abs(unadjusted_beta)-abs(adjusted_beta)'
write.csv(att,file.path(OUT,'results/cohort_level_attenuation.csv'),row.names=FALSE,quote=TRUE)
mean_rows<-list(); ii<-1
for(ref in unique(att$reference)){
 q<-att[att$reference==ref,]; n<-nrow(q); md<-mean(q$delta_abs); se<-sd(q$delta_abs)/sqrt(n); tc<-qt(.975,n-1); set.seed(20261113+ii); boot<-replicate(20000,mean(sample(q$delta_abs,n,replace=TRUE)))
 mean_rows[[ii]]<-data.frame(reference=ref,n_cohorts=n,mean_delta=md,sd_delta=sd(q$delta_abs),t_ci_low=md-tc*se,t_ci_high=md+tc*se,bootstrap_percentile_low=quantile(boot,.025),bootstrap_percentile_high=quantile(boot,.975),attenuation_n=sum(q$delta_abs>0),no_change_n=sum(q$delta_abs==0),increase_n=sum(q$delta_abs<0),bootstrap_reps=length(boot),ci_method='equal-weight cohort mean; t(7)',bootstrap_method='cohort resampling only; no inner sample bootstrap'); ii<-ii+1
}
meanatt<-do.call(rbind,mean_rows); write.csv(meanatt,file.path(OUT,'results/mean_attenuation_confidence_intervals.csv'),row.names=FALSE,quote=TRUE)
unadj<-base_meta[base_meta$effect=='unadjusted',c('reference','pooled_beta','ci_low','ci_high')]; names(unadj)[2:4]<-paste0('unadjusted_',names(unadj)[2:4]); adj<-base_meta[base_meta$effect=='adjusted',c('reference','pooled_beta','ci_low','ci_high')]; names(adj)[2:4]<-paste0('adjusted_',names(adj)[2:4]); repro<-merge(unadj,adj,by='reference'); repro$expected_unadjusted<-0.976015; repro$expected_reference_A_adjusted<-0.392659; repro$expected_reference_Bv2_adjusted<-0.604264; repro$reproduced<-ifelse(repro$reference=='Reference_A',abs(repro$adjusted_pooled_beta-repro$expected_reference_A_adjusted)<1e-6,abs(repro$adjusted_pooled_beta-repro$expected_reference_Bv2_adjusted)<1e-6)&abs(repro$unadjusted_pooled_beta-repro$expected_unadjusted)<1e-6; write.csv(repro,file.path(OUT,'results/baseline_regression_reproduction.csv'),row.names=FALSE,quote=TRUE)

fit_rows<-list(); jj<-1
for(ref in unique(fractions$reference)) for(ds in DATASETS){
 y<-fractions[fractions$reference==ref&fractions$dataset==ds,]; x<-mod[mod$dataset==ds,]; m<-merge(x,y,by='sample_id'); m$zmod<-as.numeric(scale(m$module)); m$zfrac<-as.numeric(scale(m$Myeloid_APC))
 if(nrow(m)<4 || length(unique(m$OA))<2 || length(unique(m$zfrac))<2) next
 fit<-lm(zmod~OA+zfrac,data=m); vf<-lm(OA~zfrac,data=m); vif<-1/(1-summary(vf)$r.squared); cn<-kappa(model.matrix(fit),exact=TRUE); cu<-lm(zmod~OA,data=m); ck<-cooks.distance(fit); top<-which.max(ck); g0<-m$Myeloid_APC[m$OA==0]; g1<-m$Myeloid_APC[m$OA==1]
 fit_rows[[jj]]<-data.frame(reference=ref,dataset=ds,n=nrow(m),OA_n=sum(m$OA==1),Control_n=sum(m$OA==0),composition_OA_mean=mean(g1),composition_Control_mean=mean(g0),composition_OA_median=median(g1),composition_Control_median=median(g0),composition_OA_min=min(g1),composition_OA_max=max(g1),composition_Control_min=min(g0),composition_Control_max=max(g0),composition_mean_difference=mean(g1)-mean(g0),within_cohort_spearman=suppressWarnings(cor(m$OA,m$Myeloid_APC,method='spearman')),VIF_OA=vif,condition_number=cn,unadjusted_beta=coef(cu)['OA'],unadjusted_se=coef(summary(cu))['OA','Std. Error'],adjusted_beta=coef(fit)['OA'],adjusted_se=coef(summary(fit))['OA','Std. Error'],max_cooks_distance=max(ck),max_cooks_sample=m$sample_id[top],max_cooks_OA=m$OA[top],max_cooks_Myeloid_APC=m$Myeloid_APC[top]); jj<-jj+1
}
col<-do.call(rbind,fit_rows); write.csv(col,file.path(OUT,'results/collinearity_diagnostics.csv'),row.names=FALSE,quote=TRUE)

r2<-read.csv(file.path(A3,'music_fit_by_group.csv'),stringsAsFactors=FALSE); r2<-subset(r2,run_scope=='baseline'&dataset!='ALL'&group!='Unknown'); r2sum<-aggregate(cbind(median_r2,negative_prop)~reference+dataset,r2,mean); names(r2sum)[3:4]<-c('baseline_median_r2_mean_group','baseline_negative_prop_mean_group')
low<-merge(att,r2sum,by=c('reference','dataset'),all.x=TRUE); low$low_fit_descriptive_flag<-low$baseline_median_r2_mean_group<0; low$note<-'descriptive only; no R2 threshold used for exclusion'; write.csv(low,file.path(OUT,'results/low_fit_sensitivity.csv'),row.names=FALSE,quote=TRUE)
meta_fun<-function(q){u<-rma(yi=q$unadjusted_beta,vi=q$unadjusted_se^2,method='REML',test='knha'); a<-rma(yi=q$adjusted_beta,vi=q$adjusted_se^2,method='REML',test='knha'); d<-q$delta_abs; ci<-t.test(d)$conf.int; data.frame(k=nrow(q),unadjusted_pooled_beta=as.numeric(coef(u)),unadjusted_ci_low=u$ci.lb,unadjusted_ci_high=u$ci.ub,adjusted_pooled_beta=as.numeric(coef(a)),adjusted_ci_low=a$ci.lb,adjusted_ci_high=a$ci.ub,mean_delta=mean(d),mean_delta_ci_low=ci[1],mean_delta_ci_high=ci[2],attenuation_n=sum(d>0),sign_change_n=sum(q$direction_change))}
loo_rows<-list(); kk<-1
for(ref in unique(att$reference)) for(ds in DATASETS){q<-att[att$reference==ref&att$dataset!=ds,]; z<-meta_fun(q); loo_rows[[kk]]<-data.frame(reference=ref,omitted_dataset=ds,z); kk<-kk+1}
write.csv(do.call(rbind,loo_rows),file.path(OUT,'results/leave_one_bulk_cohort_out.csv'),row.names=FALSE,quote=TRUE)

lreg<-read.csv(file.path(RES3,'LODO_adjusted_regression.csv'),stringsAsFactors=FALSE); lmeta<-read.csv(file.path(RES3,'LODO_adjusted_regression_meta.csv'),stringsAsFactors=FALSE); lreg$delta_abs<-abs(lreg$unadjusted_beta)-abs(lreg$adjusted_beta); lr<-list(); ll<-1
for(ref in unique(lreg$reference)) for(omit in unique(lreg$omit_donor[lreg$reference==ref])){q<-lreg[lreg$reference==ref&lreg$omit_donor==omit,]; d<-q$delta_abs; ti<-t.test(d)$conf.int; pm<-subset(lmeta,reference==ref&omit_donor==omit&effect=='adjusted'); lr[[ll]]<-data.frame(reference=ref,omit_donor=omit,n_cohorts=nrow(q),mean_delta=mean(d),t_ci_low=ti[1],t_ci_high=ti[2],attenuation_n=sum(d>0),no_change_n=sum(d==0),increase_n=sum(d<0),adjusted_pooled_beta=pm$pooled_beta,adjusted_ci_low=pm$ci_low,adjusted_ci_high=pm$ci_high,sign_change_n=sum(q$direction_change)); ll<-ll+1}
lodo<-do.call(rbind,lr); bl<-lodo[lodo$omit_donor=='baseline',c('reference','mean_delta')]; names(bl)[2]<-'baseline_mean_delta'; lodo<-merge(lodo,bl,by='reference'); lodo$change_from_baseline<-lodo$mean_delta-lodo$baseline_mean_delta; lodo$interpretation<-'correlated donor-dropout sensitivity runs; not independent validation'; write.csv(lodo,file.path(OUT,'results/LODO_attenuation_stability.csv'),row.names=FALSE,quote=TRUE)

pdf(file.path(OUT,'figures','Supplementary_Attenuation_Uncertainty.pdf'),width=8.5,height=5.8,useDingbats=FALSE); par(mfrow=c(1,2),mar=c(6,5,3,1),family='sans')
for(ref in unique(att$reference)){q<-att[att$reference==ref,]; plot(seq_len(nrow(q)),q$delta_abs,pch=16,col='#3B6F9E',xaxt='n',xlab='',ylab='Delta = |unadjusted beta| - |adjusted beta|',main=ref); abline(h=0,lty=2,col='grey60'); axis(1,at=seq_len(nrow(q)),labels=q$dataset,las=2,cex.axis=.65); z<-meanatt[meanatt$reference==ref,]; segments(nrow(q)+1,z$t_ci_low,nrow(q)+1,z$t_ci_high,lwd=3,col='#1F1F1F'); points(nrow(q)+1,z$mean_delta,pch=18,col='#1F1F1F'); axis(1,at=nrow(q)+1,labels='Mean',las=2,cex.axis=.7)}; dev.off()
pdf(file.path(OUT,'figures','Supplementary_Collinearity_Diagnostics.pdf'),width=8.5,height=5.8,useDingbats=FALSE); par(mfrow=c(1,2),mar=c(6,5,3,1),family='sans'); boxplot(VIF_OA~reference,data=col,col=c('#3B6F9E','#7A7A7A'),ylab='VIF for OA status',xlab='',main='Cohort-level VIF'); abline(h=5,lty=2,col='grey60'); boxplot(condition_number~reference,data=col,col=c('#3B6F9E','#7A7A7A'),ylab='Condition number',xlab='',main='Design-matrix conditioning'); dev.off()
writeLines(capture.output(sessionInfo()),file.path(OUT,'logs','sessionInfo.txt'))
writeLines(c('script=run_major_comment4_audit.R','seed=20261113','B_bootstrap=20000','primary_CI=equal-weight cohort mean with t(7)','inner_sample_bootstrap=not performed; independent-unit status incomplete','frozen_inputs=read-only reuse of revision_R1_major3 and canonical matched input'),file.path(OUT,'logs','run_parameters.txt'))
