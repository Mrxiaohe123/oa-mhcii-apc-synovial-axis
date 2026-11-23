#!/usr/bin/env Rscript
# Versioned correction package. Historical Gate directories are read-only.
# Run from the project root with: Rscript <this-script>
options(stringsAsFactors=FALSE, warn=1, scipen=999)
suppressPackageStartupMessages({
  library(metafor)
  library(Matrix)
  library(SingleCellExperiment)
  library(MuSiC)
})

ROOT <- normalizePath(getwd(), mustWork=TRUE)
BASE <- file.path(ROOT, "全部分析内容")
OUT <- file.path(BASE, "gse12021_independent_unit_correction")
RAW <- file.path(BASE, "layer25s_disease_state_recovery_and_symptom_boundary/data/raw_geo/GSE12021")
ANN <- file.path(BASE, "layer25s_disease_state_recovery_and_symptom_boundary/data/platform_annotation")
GATE1 <- file.path(BASE, "gate1_composition_linkage")
PH4 <- file.path(BASE, "phase4f_reference_based_deconvolution_sensitivity")
BV2 <- file.path(BASE, "reference_B_downstream_impact_v2")
PRIMARY <- c("GSE1919","GSE12021","GSE55235","GSE55457","GSE82107","GSE206848","GSE89408","GSE143514")
COMMON7 <- setdiff(PRIMARY, "GSE12021")
MODULE <- c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","CD74","HLA-DMA","HLA-DMB","CIITA")
PROXY <- c("LST1","TYROBP","FCER1G","AIF1","C1QA","C1QB","C1QC","CSF1R","CD68","FCGR3A","ITGAM","LYZ","SPI1","MS4A7","CTSS")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
for (d in c("source_data","scripts","logs")) dir.create(file.path(OUT,d), recursive=TRUE, showWarnings=FALSE)

LOG <- file(file.path(OUT,"logs/run.log"), open="wt")
sink(LOG, split=TRUE); sink(LOG, type="message")
cat("GSE12021 independent-unit correction; run at",format(Sys.time(),tz="Asia/Shanghai"),"\n")

sha256 <- function(f) {
  x <- system2("shasum", c("-a","256", shQuote(f)), stdout=TRUE, stderr=FALSE)
  sub("^([0-9a-fA-F]{64}).*$","\\1",x[1])
}
hash_manifest <- function() {
  targets <- c(
    file.path(BASE,"gate4_program_composition_specificity"),
    file.path(BASE,"gate4b_leave_program_out_fairness")
  )
  out <- do.call(rbind,lapply(seq_along(targets),function(i){
    d <- targets[i]
    if (!dir.exists(d)) stop("Required frozen Gate4/Gate4b directory missing: ",d)
    fs <- list.files(d,recursive=TRUE,full.names=TRUE,all.files=TRUE,no..=TRUE)
    fs <- fs[file.exists(fs) & !file.info(fs)$isdir]
    data.frame(group=if(i==1)"Gate4" else "Gate4b",path=normalizePath(fs,winslash="/"),
      size_bytes=as.numeric(file.info(fs)$size),sha256=vapply(fs,sha256,character(1)),stringsAsFactors=FALSE)
  }))
  out
}
hash_before <- hash_manifest()
write.csv(hash_before,file.path(OUT,"logs/Gate4_Gate4b_hashes_before.csv"),row.names=FALSE)

read_csv <- function(p) read.csv(p,check.names=FALSE,stringsAsFactors=FALSE,na.strings=c("","NA"))
norm_gene <- function(x) {
  x <- toupper(gsub("[[:space:]]+","",trimws(as.character(x))))
  x <- sub("///.*$|//.*$|;.*$|,.*$","",x)
  x[x %in% c("","NA","N/A","NULL","---")] <- NA_character_
  x
}
patient_from_title <- function(x) {
  m <- regexec("patient[[:space:]]+(EB[[:space:]]*[0-9]+|J[[:space:]]*[0-9]+)",x,ignore.case=TRUE)
  z <- regmatches(x,m)[[1]]
  if(length(z)<2) return(NA_character_)
  toupper(gsub("[[:space:]]+","",z[2]))
}
phenotype_from_title <- function(x) {
  if(grepl("osteoarthritis|osteoarthrisits|oateoarthritis",x,ignore.case=TRUE)) return("OA")
  if(grepl("normal",x,ignore.case=TRUE)) return("control")
  NA_character_
}

# 1) Build authoritative patient x platform map from GEO sample titles. The
# title field is used for GSM303377 because a local source_name field conflicts.
geo <- read_csv(file.path(RAW,"sample_metadata.csv"))
geo$patient_id <- vapply(geo$title,patient_from_title,character(1))
geo$phenotype <- vapply(geo$title,phenotype_from_title,character(1))
geo <- geo[!is.na(geo$patient_id) & geo$phenotype %in% c("OA","control") & geo$platform_id %in% c("GPL96","GPL97"),]
if(any(geo$patient_id=="EB140" & geo$platform_id=="GPL97" & grepl("EB149",geo$source_name_ch1))) {
  cat("Metadata conflict documented: official GEO title GSM303377=EB140; local source_name_ch1 says EB149.\n")
}
donors <- sort(unique(geo$patient_id))
if(length(donors)!=19 || sum(vapply(donors,function(id) unique(geo$phenotype[geo$patient_id==id])=="OA",logical(1)))!=10) stop("GEO donor inventory does not resolve to 10 OA + 9 controls")
map <- do.call(rbind,lapply(donors,function(id){
  d <- geo[geo$patient_id==id,]
  if(length(unique(d$phenotype))!=1) stop("Phenotype conflict for ",id)
  a <- d$geo_accession[d$platform_id=="GPL96"]; b <- d$geo_accession[d$platform_id=="GPL97"]
  if(length(a)!=1) stop("GPL96 coverage is not exactly one for ",id)
  if(length(b)>1) stop("Multiple GPL97 records for ",id)
  data.frame(patient_id=id,phenotype=unique(d$phenotype),GPL96_GSM=a,
    GPL97_GSM=if(length(b))b else "",GPL96_present=TRUE,GPL97_present=length(b)==1,
    number_of_platform_profiles=1L+length(b),
    source_metadata_evidence=paste0("NCBI GEO GSE12021 sample title(s): ",paste(d$title,collapse=" | "),
      if(id=="EB140")"; official GEO GSM303377 title resolves local source_name_ch1 EB149 discrepancy" else ""),
    stringsAsFactors=FALSE)
}))
stopifnot(sum(map$phenotype=="OA")==10,sum(map$phenotype=="control")==9,
  all(map$GPL96_present),sum(map$GPL97_present & map$phenotype=="OA")==10,
  sum(map$GPL97_present & map$phenotype=="control")==4,!anyDuplicated(map$patient_id))
write.csv(map,file.path(OUT,"GSE12021_AUTHORITATIVE_DONOR_MAP.csv"),row.names=FALSE)

# Keep a full source-record crosswalk as audit support for the correction.
raw_crosswalk <- do.call(rbind,lapply(seq_len(nrow(geo)),function(i){
  r <- geo[i,]; data.frame(GSM=r$geo_accession,patient_id=r$patient_id,phenotype=r$phenotype,
    platform=r$platform_id,canonical_GPL96=identical(r$platform_id,"GPL96"),
    canonical_analysis_include=identical(r$platform_id,"GPL96"),
    source_metadata_evidence=paste0("GEO title: ",r$title,"; source_name: ",r$source_name_ch1),stringsAsFactors=FALSE)
}))
raw_crosswalk <- raw_crosswalk[!duplicated(raw_crosswalk$GSM),]
write.csv(raw_crosswalk,file.path(OUT,"source_data/GSE12021_source_record_crosswalk.csv"),row.names=FALSE)

# 2) Frozen module gene-to-probe coverage, using precisely the historic first
# HGNC symbol rule and mean collapse across probes mapping to one gene.
coverage_for <- function(gpl) {
  a <- read_csv(file.path(ANN,paste0(gpl,".csv")))
  sym <- norm_gene(a[["Gene symbol"]])
  do.call(rbind,lapply(MODULE,function(g){
    ix <- which(sym==g)
    data.frame(platform=gpl,gene=g,present=length(ix)>0,n_probes=length(ix),
      probe_ids=if(length(ix))paste(a$ID[ix],collapse="; ") else "",
      probe_to_gene_mapping="First symbol before ///, //, semicolon or comma; uppercase and trim whitespace",
      collapsed_gene_expression_rule="Arithmetic mean of mapped probe intensities per gene, then log2(mean+1); no gene added/reweighted",
      stringsAsFactors=FALSE)
  }))
}
covA <- coverage_for("GPL96"); covB <- coverage_for("GPL97")
write.csv(covA,file.path(OUT,"GSE12021_GPL96_module_coverage.csv"),row.names=FALSE)
write.csv(transform(covB,platform_assessment=ifelse(sum(covB$present)<5,"TECHNICALLY INELIGIBLE: fixed 10-gene score coverage too sparse","eligible for secondary platform sensitivity")),
  file.path(OUT,"GSE12021_GPL97_coverage_and_sensitivity.csv"),row.names=FALSE)
if(sum(covA$present)!=10) stop("GPL96 does not permit the frozen fixed-module score; stop rather than change module")
if(sum(covB$present)>=5) stop("GPL97 eligibility threshold requires review; task predefines sparse-coverage classification")

# Helpers matching the frozen Python fixed-module effect definition.
read_probe_expr <- function(gpl_idx, map_rows) {
  f <- file.path(RAW,sprintf("series_%02d_expression_matrix.csv",gpl_idx))
  x <- read_csv(f); ids <- x[[1]]; x <- as.data.frame(lapply(x[-1],function(z)as.numeric(as.character(z))),check.names=FALSE)
  rownames(x) <- ids
  ann <- read_csv(file.path(ANN,"GPL96.csv")); sym <- norm_gene(ann[["Gene symbol"]])
  keep <- !is.na(sym) & ann$ID %in% rownames(x)
  ann <- data.frame(feature_id=ann$ID[keep],gene_symbol=sym[keep],stringsAsFactors=FALSE)
  # Match the frozen workflow: map/average probe intensities by gene, then
  # apply log2(gene_expression + 1).
  probes <- split(ann$feature_id,ann$gene_symbol)
  gx <- do.call(rbind,lapply(names(probes),function(g) colMeans(x[probes[[g]],map_rows$GPL96_GSM,drop=FALSE],na.rm=TRUE)))
  rownames(gx) <- names(probes)
  log2(gx+1)
}
z_by_gene_pop <- function(x) {
  m <- rowMeans(x,na.rm=TRUE); s <- apply(x,1,function(v)sqrt(mean((v-mean(v,na.rm=TRUE))^2,na.rm=TRUE)))
  s[!is.finite(s)|s==0] <- NA_real_
  sweep(sweep(x,1,m,"-"),1,s,"/")
}
hedges_g <- function(x,y) {
  x <- as.numeric(x[is.finite(x)]);y <- as.numeric(y[is.finite(y)])
  n1<-length(x);n2<-length(y);if(n1<2||n2<2)return(NA_real_)
  sp<-sqrt(((n1-1)*var(x)+(n2-1)*var(y))/(n1+n2-2))
  d<-(mean(x)-mean(y))/sp;j<-1-3/(4*(n1+n2)-9);d*j
}
hedges_var <- function(g,n1,n2) (n1+n2)/(n1*n2)+g^2/(2*(n1+n2-2))
boot_g_ci <- function(x,y,B=4000,seed=250704) {
  set.seed(seed); vals<-rep(NA_real_,B)
  for(i in seq_len(B)) vals[i]<-hedges_g(sample(x,length(x),replace=TRUE),sample(y,length(y),replace=TRUE))
  unname(quantile(vals,c(.025,.975),na.rm=TRUE,type=7))
}

# Exact GPL96-only primary fixed module scores (10 OA + 9 controls).
gpl96_expr <- read_probe_expr(1,map)
if(!all(MODULE%in%rownames(gpl96_expr))) stop("One or more fixed module genes missing after locked mapping")
gpl96_z <- z_by_gene_pop(gpl96_expr)
module_raw <- colMeans(gpl96_z[MODULE,,drop=FALSE],na.rm=FALSE)
module_score <- as.numeric((module_raw-mean(module_raw))/sqrt(mean((module_raw-mean(module_raw))^2)))
names(module_score)<-names(module_raw)
groups <- setNames(map$phenotype,map$GPL96_GSM)
oa <- module_score[groups[names(module_score)]=="OA"]
ctrl <- module_score[groups[names(module_score)]=="control"]
if(length(oa)!=10||length(ctrl)!=9)stop("GPL96 fixed module sample count mismatch")
g120 <- hedges_g(oa,ctrl); v120 <- hedges_var(g120,length(oa),length(ctrl)); ci120<-boot_g_ci(oa,ctrl)
effect120 <- data.frame(dataset="GSE12021",platform="GPL96",n_OA=length(oa),n_control=length(ctrl),n_total=length(oa)+length(ctrl),
  mean_OA_module_score=mean(oa),mean_control_module_score=mean(ctrl),pooled_SD=sqrt(((length(oa)-1)*var(oa)+(length(ctrl)-1)*var(ctrl))/(length(oa)+length(ctrl)-2)),
  hedges_g=g120,ci_low=ci120[1],ci_high=ci120[2],sampling_variance= v120,
  variance_formula="(n_OA+n_control)/(n_OA*n_control)+g^2/[2*(n_OA+n_control-2)]",
  module_genes_requested=length(MODULE),module_genes_present=sum(covA$present),module_scoring="within-dataset gene-wise population z-score; equal-weight mean; standardize module score across 19 profiles",
  ci_method="4000 stratified percentile bootstrap of Hedges g; seed=250704; R RNG",
  stringsAsFactors=FALSE)
write.csv(effect120,file.path(OUT,"GSE12021_GPL96_corrected_primary_effect.csv"),row.names=FALSE)

# 3) Fresh same-set Reference A/B-v2 fractions and recomputed Gate1/Gate2 rows.
load_bulk <- function(path) {
  x<-read.delim(path,check.names=FALSE,stringsAsFactors=FALSE,quote="",comment.char="")
  genes<-norm_gene(x[[1]]);mat<-as.matrix(data.frame(lapply(x[-1],function(z)as.numeric(as.character(z))),check.names=FALSE));rownames(mat)<-genes
  du<-unique(rownames(mat));mat<-do.call(rbind,lapply(du,function(g)colMeans(mat[rownames(mat)==g,,drop=FALSE],na.rm=TRUE)));rownames(mat)<-du
  mat
}
bulkA <- load_bulk(file.path(PH4,"03_deconvolution_inputs/GSE12021_bulk_expression_for_music.tsv"))
bulk_meta <- read_csv(file.path(PH4,"03_deconvolution_inputs/GSE12021_metadata_for_deconvolution.csv"))
map$bulk_sample_id <- ifelse(sub("^EB","",map$patient_id)%in%colnames(bulkA),sub("^EB","",map$patient_id),map$GPL96_GSM)
if(anyDuplicated(map$bulk_sample_id)||!all(map$bulk_sample_id%in%colnames(bulkA)))stop("Canonical 19 donor-to-expression-column mapping failed")
if(!all(map$bulk_sample_id%in%bulk_meta$sample))stop("Canonical 19 sample metadata rows missing")
map$sample_disease_match <- bulk_meta$disease_group[match(map$bulk_sample_id,bulk_meta$sample)]
if(!all(tolower(map$phenotype)==tolower(map$sample_disease_match)))stop("GEO phenotype and frozen composition metadata disagree")

run_music <- function(expr, sce, label) {
  expr[is.na(expr)]<-0
  common<-intersect(rownames(expr),rownames(sce))
  if(length(common)<1000)stop(label," has unexpectedly few common genes: ",length(common))
  fit<-MuSiC::music_prop(bulk.mtx=as.matrix(expr[common,,drop=FALSE]),sc.sce=sce[common,],clusters="coarse_class",samples="donor_id",verbose=FALSE)
  z<-as.matrix(fit$Est.prop.weighted)
  if(!setequal(rownames(z),colnames(expr))||any(!is.finite(z))||any(z<0|z>1)||max(abs(rowSums(z)-1))>1e-8)stop(label," fraction output validation failed")
  z[map$bulk_sample_id,,drop=FALSE]
}
sel <- map$bulk_sample_id
exprA <- bulkA[,sel,drop=FALSE]
sceA <- readRDS(file.path(PH4,"01_reference_construction/GSE248453_reference_primary_frozen10_excluded.rds"))
fracA <- run_music(exprA,sceA,"Reference A")

# B-v2: use the frozen, auditable B-v2 reference and exact frozen MuSiC call.
refBdir<-file.path(BV2,"inputs/reference")
genesB<-readLines(file.path(refBdir,"reference_B_genes.txt"))
mdB<-read_csv(file.path(refBdir,"reference_B_cell_metadata.csv"));rownames(mdB)<-mdB$cell_id
countsB<-Matrix::readMM(file.path(refBdir,"reference_B_counts.mtx"));rownames(countsB)<-genesB;colnames(countsB)<-mdB$cell_id
sceB<-SingleCellExperiment(assays=list(counts=countsB),colData=mdB)
bulkB<-load_bulk(file.path(BV2,"inputs/bulk/GSE12021_bulk_expression_for_music.tsv"))
if(!setequal(colnames(bulkB),colnames(bulkA))) stop("Reference A and B-v2 frozen bulk columns differ")
fracB<-run_music(bulkB[,sel,drop=FALSE],sceB,"Reference B-v2")

# Corrected module/proxy scores for Gate1/2 use same 19 GPL96 bulk profiles.
g1_module<-c(MODULE);g1_proxy<-PROXY
score_module <- function(expr,genes) {
  use<-intersect(genes,rownames(expr));if(!length(use))return(rep(NA_real_,ncol(expr)))
  z<-z_by_gene_pop(expr[use,,drop=FALSE]);colMeans(z,na.rm=TRUE)
}
modGate<-score_module(bulkA[,sel,drop=FALSE],g1_module)
proxyGate<-score_module(bulkA[,sel,drop=FALSE],g1_proxy)
names(modGate)<-sel;names(proxyGate)<-sel
oaFlag<-as.integer(map$phenotype[match(sel,map$bulk_sample_id)]=="OA")
dat<-data.frame(dataset="GSE12021",sample=sel,GSM=map$GPL96_GSM[match(sel,map$bulk_sample_id)],patient_id=map$patient_id[match(sel,map$bulk_sample_id)],
  phenotype=map$phenotype[match(sel,map$bulk_sample_id)],OA=oaFlag,module=as.numeric(modGate[sel]),proxy=as.numeric(proxyGate[sel]),
  Myeloid_APC_A=fracA[,"Myeloid_APC"],Myeloid_APC_Bv2=fracB[,"Myeloid_APC"],stringsAsFactors=FALSE)
if(nrow(dat)!=19||anyDuplicated(dat$patient_id)||sum(dat$OA)!=10||sum(dat$OA==0)!=9)stop("Corrected Gate1/2 sample set must be 19 unique participants (10/9)")

lm_row <- function(formula,d,term="OA") {
  fit<-lm(formula,data=d);cc<-summary(fit)$coefficients[term,];df<-fit$df.residual
  c(beta=unname(cc["Estimate"]),se=unname(cc["Std. Error"]),low=unname(cc["Estimate"]-qt(.975,df)*cc["Std. Error"]),
    high=unname(cc["Estimate"]+qt(.975,df)*cc["Std. Error"]),p=unname(cc["Pr(>|t|)"]),df=df)
}
g1_metrics <- function(d,ref) {
  y<-d[[paste0("Myeloid_APC_",ref)]];rho<-cor(d$module,y,method="spearman")
  oa_y<-y[d$OA==1];co_y<-y[d$OA==0]
  es<-metafor::escalc(measure="SMD",m1i=mean(oa_y),sd1i=sd(oa_y),n1i=length(oa_y),m2i=mean(co_y),sd2i=sd(co_y),n2i=length(co_y),append=TRUE)
  u<-lm_row(scale(module)~OA,d);a<-lm_row(scale(module)~OA+scale(y),d)
  proxy_adj<-lm_row(scale(module)~OA+scale(proxy),d)
  data.frame(dataset="GSE12021",reference=ref,n=nrow(d),n_OA=sum(d$OA),n_control=sum(d$OA==0),
    spearman_rho=rho,Fisher_z=atanh(rho),Fisher_variance=1/(nrow(d)-3),
    module_genes_available=length(intersect(MODULE,rownames(bulkA))),
    oa_myeloid_mean_OA=mean(oa_y),oa_myeloid_mean_control=mean(co_y),oa_myeloid_difference=mean(oa_y)-mean(co_y),
    OA_APC_hedges_g=as.numeric(es$yi),OA_APC_variance=as.numeric(es$vi),OA_APC_SE=sqrt(as.numeric(es$vi)),
    unadjusted_beta=u["beta"],unadjusted_SE=u["se"],unadjusted_low=u["low"],unadjusted_high=u["high"],
    adjusted_beta=a["beta"],adjusted_SE=a["se"],adjusted_low=a["low"],adjusted_high=a["high"],adjusted_df=a["df"],
    attenuation=abs(a["beta"])<abs(u["beta"]),attenuation_class=ifelse(abs(a["beta"])<abs(u["beta"]),"attenuated",ifelse(abs(a["beta"])>abs(u["beta"]),"amplified","unchanged")),
    direction_change=sign(a["beta"])!=sign(u["beta"]),marker_proxy_genes_available=length(intersect(PROXY,rownames(bulkA))),
    marker_proxy_adjusted_beta=proxy_adj["beta"],marker_proxy_adjusted_SE=proxy_adj["se"],marker_proxy_adjusted_low=proxy_adj["low"],marker_proxy_adjusted_high=proxy_adj["high"],
    marker_proxy_attenuation=abs(proxy_adj["beta"])<abs(u["beta"]),stringsAsFactors=FALSE)
}
g1row<-g1_metrics(dat,"A");g2row<-g1_metrics(dat,"Bv2")
if(any(abs(c(g1row$spearman_rho,g2row$spearman_rho))>1))stop("Correlation bounds violated")
write.csv(dat,file.path(OUT,"GSE12021_GPL96_corrected_scores.csv"),row.names=FALSE)
write.csv(g1row,file.path(OUT,"gate1_GSE12021_corrected.csv"),row.names=FALSE)
write.csv(fracB,file.path(OUT,"Bv2_GSE12021_corrected_fractions.csv"),row.names=TRUE)
write.csv(g2row,file.path(OUT,"gate2_Bv2_GSE12021_corrected.csv"),row.names=FALSE)

# 4) REML/HK primary meta-analysis; seven historical dataset rows must stay
# exactly as frozen, with only the GSE12021 effect/n/variance replaced.
old_primary<-read_csv(file.path(BASE,"outputs_reml_hk_robustness/study_level_frozen_inputs.csv"))
source_meta<-read_csv(file.path(BASE,"layer25s_disease_state_recovery_and_symptom_boundary/01_legacy_meta/study_level_module_effects.csv"))
old8<-source_meta[source_meta$dataset_id%in%PRIMARY & source_meta$module=="MHCII_antigen_presentation_module",]
old8<-old8[match(PRIMARY,old8$dataset_id),]
if(nrow(old8)!=8)stop("Frozen eight study rows unavailable")
old7<-old8[old8$dataset_id%in%COMMON7,]
old7_orig<-old_primary[old_primary$dataset_id%in%COMMON7,]
old7_orig<-old7_orig[match(COMMON7,old7_orig$dataset_id),]
old7<-old7[match(COMMON7,old7$dataset_id),]
other7_check<-data.frame(dataset=COMMON7,
  old_effect=old7$hedges_g,correction_effect=old7_orig$hedges_g,
  old_n_OA=old7$n_oa,correction_n_OA=old7_orig$n_oa,
  old_n_control=old7$n_control,correction_n_control=old7_orig$n_control,
  old_variance=old7$variance,correction_variance=old7_orig$variance,
  unchanged=abs(old7$hedges_g-old7_orig$hedges_g)<1e-12 & old7$n_oa==old7_orig$n_oa & old7$n_control==old7_orig$n_control & abs(old7$variance-old7_orig$variance)<1e-12,
  stringsAsFactors=FALSE)
write.csv(other7_check,file.path(OUT,"other_seven_primary_input_validation.csv"),row.names=FALSE)
if(!all(other7_check$unchanged)) stop("Other-seven primary study inputs changed")
preplot_hashes<-read_csv(file.path(ROOT,"manuscript_reconstruction/preplot_reconciliation/PREPLOT_QA_MANIFEST.csv"))
source_row<-preplot_hashes[basename(preplot_hashes$path)=="study_level_module_effects.csv",]
source_file<-file.path(BASE,"layer25s_disease_state_recovery_and_symptom_boundary/01_legacy_meta/study_level_module_effects.csv")
source_hash_current<-sha256(source_file)
if(nrow(source_row)!=1 || source_hash_current!=source_row$sha256) stop("Frozen source hash differs from pre-correction QA manifest")
write.csv(data.frame(path=normalizePath(source_file,winslash="/"),expected_sha256=source_row$sha256,current_sha256=source_hash_current,unchanged=source_hash_current==source_row$sha256),
  file.path(OUT,"logs/other_seven_primary_source_hash_validation.csv"),row.names=FALSE)
new8<-old8
new8$n_oa[new8$dataset_id=="GSE12021"]<-10L
new8$n_control[new8$dataset_id=="GSE12021"]<-9L
new8$hedges_g[new8$dataset_id=="GSE12021"]<-g120
new8$variance[new8$dataset_id=="GSE12021"]<-v120
new8$sampling_variance_recomputed<-with(new8,(n_oa+n_control)/(n_oa*n_control)+hedges_g^2/(2*(n_oa+n_control-2)))
if(any(abs(new8$variance-new8$sampling_variance_recomputed)>1e-12))stop("Primary variance formula mismatch")
if(!all(vapply(COMMON7,function(ds)identical(as.numeric(new8$hedges_g[new8$dataset_id==ds]),as.numeric(old8$hedges_g[old8$dataset_id==ds])),logical(1))))stop("One of the other seven effects changed")

pool_hk <- function(yi,vi,label,correlation=FALSE) {
  ok<-is.finite(yi)&is.finite(vi)&vi>0;yi<-yi[ok];vi<-vi[ok];k<-length(yi)
  if(k<3)stop("REML/HK pooling requires >=3 datasets: ",label)
  fit<-rma(yi=yi,vi=vi,method="REML",test="knha",control=list(maxiter=1000))
  pr<-predict(fit);vals<-c(as.numeric(fit$b),fit$ci.lb,fit$ci.ub,pr$pi.lb,pr$pi.ub)
  if(correlation) vals<-tanh(vals)
  if(correlation && (any(!is.finite(vals))||any(vals< -1|vals>1)))stop("Fisher-z back-transform bounds failed: ",label)
  data.frame(metric=label,k=k,pooled_effect=vals[1],ci_low=vals[2],ci_high=vals[3],pi_low=vals[4],pi_high=vals[5],
    p_value=as.numeric(fit$pval),HK_SE=as.numeric(fit$se),HK_t=as.numeric(fit$zval),df=as.numeric(fit$ddf),tau2=as.numeric(fit$tau2),I2=as.numeric(fit$I2),Q=as.numeric(fit$QE),
    inference="REML + Hartung-Knapp; t(k-1)",prediction_inference="Knapp-Hartung t(k-1)",scale=if(correlation)"Fisher-z fit; all displayed effects and intervals tanh-back-transformed" else "Hedges g",stringsAsFactors=FALSE)
}
pool_primary <- function(d,label="primary_REML_HK") pool_hk(d$hedges_g,d$sampling_variance_recomputed,label)
primary_sum<-pool_primary(new8)
write.csv(primary_sum,file.path(OUT,"primary_meta_GSE12021_corrected.csv"),row.names=FALSE)
new7<-new8[new8$dataset_id%in%COMMON7,]
without12021<-pool_primary(new7,"primary_without_GSE12021_REML_HK")
without12021$positive_datasets<-sum(new7$hedges_g>0);without12021$total_datasets<-nrow(new7);without12021$direction_consistency<-paste0(without12021$positive_datasets,"/",without12021$total_datasets," positive")
write.csv(without12021,file.path(OUT,"primary_meta_without_GSE12021.csv"),row.names=FALSE)
loo_primary<-do.call(rbind,lapply(PRIMARY,function(drop){d<-new8[new8$dataset_id!=drop,];q<-pool_primary(d,paste0("omit_",drop));q$omitted_dataset<-drop;q$pooled_direction_positive<-q$pooled_effect>0;q$hk_ci_crosses_zero<-q$ci_low<=0&&q$ci_high>=0;q}))
write.csv(loo_primary,file.path(OUT,"primary_meta_LODO_GSE12021_corrected.csv"),row.names=FALSE)
write.csv(new8,file.path(OUT,"primary_meta_study_inputs_corrected.csv"),row.names=FALSE)

# 5) Corrected Gate1 & B-v2 pooled summaries. Reuse all seven frozen dataset
# rows and replace GSE12021; Fisher-z is the fitted correlation scale.
oldA_direct<-read_csv(file.path(GATE1,"module_vs_myeloid_dataset_results.csv"));oldA_direct<-oldA_direct[oldA_direct$role=="primary"&oldA_direct$dataset%in%COMMON7,]
oldA_comp<-read_csv(file.path(GATE1,"oa_vs_control_myeloid_results.csv"));oldA_comp<-oldA_comp[oldA_comp$role=="primary"&oldA_comp$dataset%in%COMMON7,]
oldA_att<-read_csv(file.path(GATE1,"oa_attenuation_dataset_results.csv"));oldA_att<-oldA_att[oldA_att$role=="primary"&oldA_att$dataset%in%COMMON7,]
oldA_mark<-read_csv(file.path(GATE1,"marker_proxy_sensitivity_results.csv"));oldA_mark<-oldA_mark[oldA_mark$role=="primary"&oldA_mark$dataset%in%COMMON7,]
oldB_direct<-read_csv(file.path(BV2,"gate2_Bv2_module_composition.csv"));oldB_direct<-oldB_direct[oldB_direct$reference_version=="Bv2"&oldB_direct$role=="primary"&oldB_direct$dataset%in%COMMON7,]
oldB_comp<-read_csv(file.path(BV2,"gate2_Bv2_oa_composition.csv"));oldB_comp<-oldB_comp[oldB_comp$reference_version=="Bv2"&oldB_comp$role=="primary"&oldB_comp$dataset%in%COMMON7,]
oldB_att<-read_csv(file.path(BV2,"gate2_Bv2_attenuation.csv"));oldB_att<-oldB_att[oldB_att$reference_version=="Bv2"&oldB_att$role=="primary"&oldB_att$dataset%in%COMMON7,]
if(any(lengths(list(oldA_direct$dataset,oldA_comp$dataset,oldA_att$dataset,oldA_mark$dataset,oldB_direct$dataset,oldB_comp$dataset,oldB_att$dataset))!=7))stop("Frozen seven-dataset Gate1/B-v2 rows incomplete")

pool_gate1 <- function() {
  rows<-list()
  allD<-rbind(oldA_direct[,c("dataset","n","spearman_rho")],data.frame(dataset="GSE12021",n=19,spearman_rho=g1row$spearman_rho))
  rows[[length(rows)+1]]<-pool_hk(atanh(allD$spearman_rho),1/(allD$n-3),"Reference_A_module_composition_rho",TRUE)
  allC<-rbind(data.frame(dataset=oldA_comp$dataset,yi=oldA_comp$hedges_g,vi=oldA_comp$hedges_se^2),data.frame(dataset="GSE12021",yi=g1row$OA_APC_hedges_g,vi=g1row$OA_APC_variance))
  rows[[length(rows)+1]]<-pool_hk(allC$yi,allC$vi,"Reference_A_OA_APC_Hedges_g")
  for(m in c("unadjusted","adjusted")){
    b<-if(m=="unadjusted")"model1_beta" else "model3_beta";lo<-if(m=="unadjusted")"model1_low" else "model3_low";hi<-if(m=="unadjusted")"model1_high" else "model3_high";df<-if(m=="unadjusted")oldA_att$n-2 else oldA_att$n-3
    old_se<-(oldA_att[[hi]]-oldA_att[[lo]])/(2*qt(.975,df))
    newb<-if(m=="unadjusted")g1row$unadjusted_beta else g1row$adjusted_beta;newse<-if(m=="unadjusted")g1row$unadjusted_SE else g1row$adjusted_SE
    vals<-c(oldA_att[[b]],newb);ses<-c(old_se,newse)
    rows[[length(rows)+1]]<-pool_hk(vals,ses^2,paste0("Reference_A_OA_module_beta_",m))
  }
  for(m in c("unadjusted","adjusted")){
    b<-if(m=="unadjusted")"model1_beta" else "model3_beta";lo<-if(m=="unadjusted")"model1_low" else "model3_low";hi<-if(m=="unadjusted")"model1_high" else "model3_high";df<-if(m=="unadjusted")oldA_mark$n-2 else oldA_mark$n-3
    old_se<-(oldA_mark[[hi]]-oldA_mark[[lo]])/(2*qt(.975,df))
    newb<-if(m=="unadjusted")g1row$unadjusted_beta else g1row$marker_proxy_adjusted_beta
    newse<-if(m=="unadjusted")g1row$unadjusted_SE else g1row$marker_proxy_adjusted_SE
    rows[[length(rows)+1]]<-pool_hk(c(oldA_mark[[b]],newb),c(old_se,newse)^2,paste0("Reference_A_marker_proxy_beta_",m))
  }
  do.call(rbind,rows)
}
gate1_pool<-pool_gate1();write.csv(gate1_pool,file.path(OUT,"gate1_corrected_pooled_summary.csv"),row.names=FALSE)

pool_bv2 <- function() {
  rows<-list();allD<-rbind(oldB_direct[,c("dataset","n","spearman_rho")],data.frame(dataset="GSE12021",n=19,spearman_rho=g2row$spearman_rho))
  rows[[1]]<-pool_hk(atanh(allD$spearman_rho),1/(allD$n-3),"Bv2_module_composition_rho",TRUE)
  allC<-rbind(data.frame(dataset=oldB_comp$dataset,yi=oldB_comp$hedges_g,vi=oldB_comp$variance),data.frame(dataset="GSE12021",yi=g2row$OA_APC_hedges_g,vi=g2row$OA_APC_variance))
  rows[[2]]<-pool_hk(allC$yi,allC$vi,"Bv2_OA_APC_Hedges_g")
  for(m in c("unadjusted","adjusted")){
    b<-if(m=="unadjusted")"unadjusted_beta" else "adjusted_beta";se<-if(m=="unadjusted")"unadjusted_SE" else "adjusted_SE"
    nb<-if(m=="unadjusted")g2row$unadjusted_beta else g2row$adjusted_beta;ns<-if(m=="unadjusted")g2row$unadjusted_SE else g2row$adjusted_SE
    rows[[length(rows)+1]]<-pool_hk(c(oldB_att[[b]],nb),c(oldB_att[[se]],ns)^2,paste0("Bv2_OA_module_beta_",m))
  }
  do.call(rbind,rows)
}
gate2_pool<-pool_bv2();write.csv(gate2_pool,file.path(OUT,"gate2_Bv2_corrected_pooled_summary.csv"),row.names=FALSE)

# Cross-reference directions and attenuation classification per dataset.
ad<-read_csv(file.path(GATE1,"module_vs_myeloid_dataset_results.csv"));ac<-read_csv(file.path(GATE1,"oa_vs_control_myeloid_results.csv"));aa<-read_csv(file.path(GATE1,"oa_attenuation_dataset_results.csv"))
bd<-read_csv(file.path(BV2,"gate2_Bv2_module_composition.csv"));bc<-read_csv(file.path(BV2,"gate2_Bv2_oa_composition.csv"));ba<-read_csv(file.path(BV2,"gate2_Bv2_attenuation.csv"))
conc<-do.call(rbind,lapply(PRIMARY,function(ds){
  if(ds=="GSE12021") {rhoA<-g1row$spearman_rho;rhoB<-g2row$spearman_rho;diffA<-g1row$oa_myeloid_difference;diffB<-g2row$oa_myeloid_difference;attenA<-g1row$attenuation;attenB<-g2row$attenuation}
  else {
    x<-ad[ad$dataset==ds,];y<-bd[bd$dataset==ds&bd$reference_version=="Bv2",];u<-ac[ac$dataset==ds,];v<-bc[bc$dataset==ds&bc$reference_version=="Bv2",];q<-aa[aa$dataset==ds,];r<-ba[ba$dataset==ds&ba$reference_version=="Bv2",]
    rhoA<-x$spearman_rho;rhoB<-y$spearman_rho;diffA<-u$raw_difference;diffB<-v$raw_difference;attenA<-q$classification=="attenuated";attenB<-r$attenuates
  }
  data.frame(dataset=ds,Reference_A_rho=rhoA,Reference_A_rho_direction=ifelse(rhoA>0,"positive","negative"),Bv2_rho=rhoB,Bv2_rho_direction=ifelse(rhoB>0,"positive","negative"),
    Reference_A_OA_composition_difference=diffA,Reference_A_OA_composition_direction=ifelse(diffA>0,"higher_OA","lower_OA"),
    Bv2_OA_composition_difference=diffB,Bv2_OA_composition_direction=ifelse(diffB>0,"higher_OA","lower_OA"),
    Reference_A_attenuation=ifelse(attenA,"attenuated","not_attenuated"),Bv2_attenuation=ifelse(attenB,"attenuated","not_attenuated"),
    module_rho_concordant=(rhoA>0)==(rhoB>0),OA_composition_concordant=(diffA>0)==(diffB>0),attenuation_concordant=attenA==attenB,
    overall_concordance=((rhoA>0)==(rhoB>0))&&((diffA>0)==(diffB>0))&&(attenA==attenB),stringsAsFactors=FALSE)
}))
write.csv(conc,file.path(OUT,"reference_A_Bv2_corrected_concordance.csv"),row.names=FALSE)

# 6) Old-vs-corrected results and corrected sample accounting.
oldPrimarySummary<-read_csv(file.path(BASE,"outputs_reml_hk_robustness/reml_hk_primary_summary.csv"))
oldA<-read_csv(file.path(GATE1,"pooled_meta_results.csv"));oldB<-read_csv(file.path(BV2,"gate2_Bv2_pooled_summary.csv"));oldB<-oldB[oldB$reference_version=="Bv2",]
oldA_verified_all<-read_csv(file.path(BASE,"gate4_program_composition_specificity/logs/index_meta_reproduction.csv"))
oldA_verified<-oldA_verified_all[oldA_verified_all$reference=="A"&oldA_verified_all$analysis=="legacy_index_eight_dataset_diagnostic",]
if(nrow(oldA_verified)!=1)stop("Previously audited Reference A Fisher-z result is unavailable")
old_key<-list(
 data.frame(domain="Primary OA meta-analysis",metric="pooled_Hedges_g",old_estimate=oldPrimarySummary$pooled_effect,new_estimate=primary_sum$pooled_effect,old_ci_low=oldPrimarySummary$ci_low,old_ci_high=oldPrimarySummary$ci_high,new_ci_low=primary_sum$ci_low,new_ci_high=primary_sum$ci_high,old_pi_low=oldPrimarySummary$prediction_low,old_pi_high=oldPrimarySummary$prediction_high,new_pi_low=primary_sum$pi_low,new_pi_high=primary_sum$pi_high,old_I2=oldPrimarySummary$I2_reml,new_I2=primary_sum$I2,stringsAsFactors=FALSE),
 data.frame(domain="GSE12021 primary",metric="Hedges_g",old_estimate=old8$hedges_g[old8$dataset_id=="GSE12021"],new_estimate=g120,old_ci_low=old8$ci_low[old8$dataset_id=="GSE12021"],old_ci_high=old8$ci_high[old8$dataset_id=="GSE12021"],new_ci_low=ci120[1],new_ci_high=ci120[2],old_pi_low=NA,old_pi_high=NA,new_pi_low=NA,new_pi_high=NA,old_I2=NA,new_I2=NA),
 data.frame(domain="Reference A",metric="module_composition_rho",old_estimate=oldA_verified$pooled,new_estimate=gate1_pool$pooled_effect[gate1_pool$metric=="Reference_A_module_composition_rho"],old_ci_low=oldA_verified$ci_low,old_ci_high=oldA_verified$ci_high,new_ci_low=gate1_pool$ci_low[gate1_pool$metric=="Reference_A_module_composition_rho"],new_ci_high=gate1_pool$ci_high[gate1_pool$metric=="Reference_A_module_composition_rho"],old_pi_low=oldA_verified$pi_low,old_pi_high=oldA_verified$pi_high,new_pi_low=gate1_pool$pi_low[gate1_pool$metric=="Reference_A_module_composition_rho"],new_pi_high=gate1_pool$pi_high[gate1_pool$metric=="Reference_A_module_composition_rho"],old_I2=oldA_verified$I2,new_I2=gate1_pool$I2[gate1_pool$metric=="Reference_A_module_composition_rho"]),
 data.frame(domain="B-v2",metric="module_composition_rho",old_estimate=oldB$pooled[oldB$metric=="module_composition_rho"],new_estimate=gate2_pool$pooled_effect[gate2_pool$metric=="Bv2_module_composition_rho"],old_ci_low=oldB$ci_low[oldB$metric=="module_composition_rho"],old_ci_high=oldB$ci_high[oldB$metric=="module_composition_rho"],new_ci_low=gate2_pool$ci_low[gate2_pool$metric=="Bv2_module_composition_rho"],new_ci_high=gate2_pool$ci_high[gate2_pool$metric=="Bv2_module_composition_rho"],old_pi_low=oldB$pi_low[oldB$metric=="module_composition_rho"],old_pi_high=oldB$pi_high[oldB$metric=="module_composition_rho"],new_pi_low=gate2_pool$pi_low[gate2_pool$metric=="Bv2_module_composition_rho"],new_pi_high=gate2_pool$pi_high[gate2_pool$metric=="Bv2_module_composition_rho"],old_I2=oldB$I2[oldB$metric=="module_composition_rho"],new_I2=gate2_pool$I2[gate2_pool$metric=="Bv2_module_composition_rho"]),
 data.frame(domain="Reference A attenuation",metric="attenuating_datasets",old_estimate=NA,new_estimate=sum(aa$dataset%in%COMMON7 & aa$role=="primary" & aa$classification=="attenuated")+as.integer(g1row$attenuation),old_ci_low=NA,old_ci_high=NA,new_ci_low=NA,new_ci_high=NA,old_pi_low=NA,old_pi_high=NA,new_pi_low=NA,new_pi_high=NA,old_I2=NA,new_I2=NA),
 data.frame(domain="B-v2 attenuation",metric="attenuating_datasets",old_estimate=NA,new_estimate=sum(ba$dataset%in%COMMON7 & ba$role=="primary" & ba$reference_version=="Bv2" & ba$attenuates)+as.integer(g2row$attenuation),old_ci_low=NA,old_ci_high=NA,new_ci_low=NA,new_ci_high=NA,old_pi_low=NA,old_pi_high=NA,new_pi_low=NA,new_pi_high=NA,old_I2=NA,new_I2=NA)
)
oldvs<-do.call(rbind,old_key)
oldA_att_all<-read_csv(file.path(GATE1,"oa_attenuation_dataset_results.csv"))
oldB_att_all<-read_csv(file.path(BV2,"gate2_Bv2_attenuation.csv"))
oldvs$old_estimate[oldvs$domain=="Reference A attenuation"]<-sum(oldA_att_all$role=="primary"&oldA_att_all$classification=="attenuated")
oldvs$old_estimate[oldvs$domain=="B-v2 attenuation"]<-sum(oldB_att_all$role=="primary"&oldB_att_all$reference_version=="Bv2"&oldB_att_all$attenuates)
write.csv(oldvs,file.path(OUT,"old_vs_corrected_key_results.csv"),row.names=FALSE)

old_meta_total<-sum(old8$n_oa+old8$n_control);new_meta_total<-sum(new8$n_oa+new8$n_control)
other <- new8[new8$dataset_id%in%COMMON7,]
gate1rows<-read_csv(file.path(GATE1,"source_data/oa_attenuation_dataset_results.csv")); gate1rows<-gate1rows[gate1rows$role=="primary"&gate1rows$dataset%in%COMMON7,]
account <- data.frame(dataset=PRIMARY,OA_analysis_units=new8$n_oa[match(PRIMARY,new8$dataset_id)],control_analysis_units=new8$n_control[match(PRIMARY,new8$dataset_id)],
  total_analysis_units=new8$n_oa[match(PRIMARY,new8$dataset_id)]+new8$n_control[match(PRIMARY,new8$dataset_id)],
  biological_donors_verified=ifelse(PRIMARY=="GSE12021",19,NA),unit_label=ifelse(PRIMARY=="GSE12021","GPL96 unique participant profiles; one per GEO patient ID","Frozen study-level expression-analysis units; donor uniqueness not independently verified"),stringsAsFactors=FALSE)
account$old_total_profile_records<-c(old8$n_oa+old8$n_control)[match(PRIMARY,old8$dataset_id)]
account$removed_duplicate_platform_records<-account$old_total_profile_records-account$total_analysis_units
account$common_seven_comparator_units<-ifelse(account$dataset%in%COMMON7,account$total_analysis_units,0L)
write.csv(account,file.path(OUT,"corrected_primary_dataset_accounting.csv"),row.names=FALSE)
write.csv(account,file.path(ROOT,"manuscript_reconstruction/PRIMARY_DATASET_ACCOUNTING.csv"),row.names=FALSE)

# Old vs corrected sample count total audit.
cat("Primary analysis units old:",old_meta_total," corrected:",new_meta_total," removed:",old_meta_total-new_meta_total,"\n")

# Replay all eight old Gate1/B-v2 direct-correlation values and pooled estimates.
oldA_replay_rows<-ad[ad$role=="primary"&ad$dataset%in%PRIMARY,]
oldB_replay_rows<-bd[bd$reference_version=="Bv2"&bd$role=="primary"&bd$dataset%in%PRIMARY,]
if(nrow(oldA_replay_rows)!=8||nrow(oldB_replay_rows)!=8)stop("Old correlation replay requires the full frozen eight-dataset set")
oldA_replay<-pool_hk(atanh(oldA_replay_rows$spearman_rho),1/(oldA_replay_rows$n-3),"primary_module_vs_myeloid_spearman",TRUE)
oldA_frozen<-oldA[oldA$analysis=="primary_module_vs_myeloid_spearman",]
oldB_replay<-pool_hk(atanh(oldB_replay_rows$spearman_rho),1/(oldB_replay_rows$n-3),"module_composition_rho",TRUE)
oldB_frozen<-read_csv(file.path(BV2,"gate2_Bv2_pooled_summary.csv"));oldB_frozen<-oldB_frozen[oldB_frozen$reference_version=="Bv2"&oldB_frozen$metric=="module_composition_rho",]
replay_checks<-data.frame(analysis=c("Gate1 old A rho: independent Fisher-z/REML-HK replay","B-v2 old rho pooled replay"),
  independently_replayed=c(oldA_replay$pooled_effect,oldB_replay$pooled_effect),
  frozen_table_value=c(oldA_frozen$pooled_effect,oldB_frozen$pooled),
  difference_from_frozen_table=abs(c(oldA_replay$pooled_effect,oldB_replay$pooled_effect)-c(oldA_frozen$pooled_effect,oldB_frozen$pooled)),
  verified_reference=c(0.514375,oldB_frozen$pooled),
  tolerance=c(1e-6,1e-8),
  pass=c(abs(oldA_replay$pooled_effect-0.514375)<1e-6,abs(oldB_replay$pooled_effect-oldB_frozen$pooled)<1e-8),
  interpretation=c("The historical Gate1 pooled_meta_results.csv is stale/non-reproducing; dataset-level frozen inputs independently reproduce the already audited Fisher-z result.","B-v2 dataset-level values reproduce the frozen pooled summary."),stringsAsFactors=FALSE)
write.csv(replay_checks,file.path(OUT,"logs/meta_replay_checks.csv"),row.names=FALSE)
if(!all(replay_checks$pass))stop("At least one independent Fisher-z meta-analysis validation failed")

# 7) Scientific impact classification; thresholds/rationale are recorded in README.
oldA_pos<-sum(ad$role=="primary"&ad$spearman_rho>0);newA_pos<-sum(conc$Reference_A_rho>0)
oldBsum<-read_csv(file.path(BV2,"gate2_Bv2_pooled_summary.csv"));oldBsum<-oldBsum[oldBsum$reference_version=="Bv2",]
oldB_pos<-sum(bd$role=="primary"&bd$reference_version=="Bv2"&bd$spearman_rho>0);newB_pos<-sum(conc$Bv2_rho>0)
oldA_att_n<-sum(oldA_att_all$role=="primary"&oldA_att_all$classification=="attenuated")
oldB_att_n<-sum(oldB_att_all$role=="primary"&oldB_att_all$reference_version=="Bv2"&oldB_att_all$attenuates)
newA_att_n<-sum(conc$Reference_A_attenuation=="attenuated");newB_att_n<-sum(conc$Bv2_attenuation=="attenuated")
primary_ci_positive<-primary_sum$ci_low>0
impact_primary<-if(primary_sum$pooled_effect<=0)"LEVEL D — CENTRAL CONCLUSION CHANGES" else if(!primary_ci_positive)"LEVEL C — MATERIAL PARTIAL IMPACT" else if(abs(primary_sum$pooled_effect-oldPrimarySummary$pooled_effect)>.1||abs(primary_sum$I2-oldPrimarySummary$I2_reml)>5)"LEVEL B — QUANTITATIVE CHANGE, SAME CONCLUSIONS" else "LEVEL A — NO MATERIAL IMPACT"
impact_A<-if(gate1_pool$pooled_effect[gate1_pool$metric=="Reference_A_module_composition_rho"]<=0)"LEVEL D — CENTRAL CONCLUSION CHANGES" else if(gate1_pool$ci_low[gate1_pool$metric=="Reference_A_module_composition_rho"]<=0)"LEVEL C — MATERIAL PARTIAL IMPACT" else if(abs(gate1_pool$pooled_effect[gate1_pool$metric=="Reference_A_module_composition_rho"]-oldA_verified$pooled)>.1||oldA_pos!=newA_pos)"LEVEL B — QUANTITATIVE CHANGE, SAME CONCLUSIONS" else "LEVEL A — NO MATERIAL IMPACT"
impact_B<-if(gate2_pool$pooled_effect[gate2_pool$metric=="Bv2_module_composition_rho"]<=0)"LEVEL D — CENTRAL CONCLUSION CHANGES" else if(gate2_pool$ci_low[gate2_pool$metric=="Bv2_module_composition_rho"]<=0)"LEVEL C — MATERIAL PARTIAL IMPACT" else if(abs(gate2_pool$pooled_effect[gate2_pool$metric=="Bv2_module_composition_rho"]-oldBsum$pooled[oldBsum$metric=="module_composition_rho"])>.1||oldB_pos!=newB_pos)"LEVEL B — QUANTITATIVE CHANGE, SAME CONCLUSIONS" else "LEVEL A — NO MATERIAL IMPACT"
impact_atten<-if(newA_att_n==0&&newB_att_n==0)"LEVEL C — MATERIAL PARTIAL IMPACT" else if(newA_att_n!=oldA_att_n||newB_att_n!=oldB_att_n)"LEVEL B — QUANTITATIVE CHANGE, SAME CONCLUSIONS" else "LEVEL A — NO MATERIAL IMPACT"
overall<-if(grepl("LEVEL D",impact_primary)||grepl("LEVEL D",impact_A)||grepl("LEVEL D",impact_B))"LEVEL D — CENTRAL CONCLUSION CHANGES" else if(grepl("LEVEL C",impact_primary)||grepl("LEVEL C",impact_A)||grepl("LEVEL C",impact_B))"LEVEL C — MATERIAL PARTIAL IMPACT" else if(any(grepl("LEVEL B",c(impact_primary,impact_A,impact_B,impact_atten))))"LEVEL B — QUANTITATIVE CHANGE, SAME CONCLUSIONS" else "LEVEL A — NO MATERIAL IMPACT"
impact<-data.frame(domain=c("Primary OA meta-analysis","Reference A composition linkage","B-v2 composition linkage","Composition attenuation","Overall manuscript story"),
  classification=c(impact_primary,impact_A,impact_B,impact_atten,overall),
  basis=c(paste0("Corrected pooled g=",signif(primary_sum$pooled_effect,5),"; 95% CI ",signif(primary_sum$ci_low,5)," to ",signif(primary_sum$ci_high,5),"; seven-dataset omit estimate=",signif(without12021$pooled_effect,5)),
    paste0("A rho pooled=",signif(gate1_pool$pooled_effect[gate1_pool$metric=="Reference_A_module_composition_rho"],5),"; positive datasets ",oldA_pos,"→",newA_pos),
    paste0("B-v2 rho pooled=",signif(gate2_pool$pooled_effect[gate2_pool$metric=="Bv2_module_composition_rho"],5),"; positive datasets ",oldB_pos,"→",newB_pos),
    paste0("attenuation A ",oldA_att_n,"→",newA_att_n,"/8; B-v2 ",oldB_att_n,"→",newB_att_n,"/8"),
    paste0("Composite classification follows whether corrected primary and both composition links preserve positive pooled direction and interval interpretation; Gate4/4b remain frozen and unaffected.")),
  stringsAsFactors=FALSE)
write.csv(impact,file.path(OUT,"scientific_impact_assessment.csv"),row.names=FALSE)

# 8) Hard validation and before/after checksums.
stopifnot(nrow(map)==19,sum(map$phenotype=="OA")==10,sum(map$phenotype=="control")==9,
  !anyDuplicated(map$patient_id),all(map$GPL96_present),length(intersect(map$GPL97_GSM,map$GPL96_GSM))==0,
  setequal(dat$GSM,map$GPL96_GSM),setequal(dat$sample,sel),sum(new8$n_oa+new8$n_control)==158)
hash_after<-hash_manifest();write.csv(hash_after,file.path(OUT,"logs/Gate4_Gate4b_hashes_after.csv"),row.names=FALSE)
mergehash<-merge(hash_before,hash_after,by=c("group","path","size_bytes"),suffixes=c("_before","_after"),all=TRUE)
mergehash$unchanged<-!is.na(mergehash$sha256_before)&!is.na(mergehash$sha256_after)&mergehash$sha256_before==mergehash$sha256_after
write.csv(mergehash,file.path(OUT,"logs/Gate4_Gate4b_checksum_validation.csv"),row.names=FALSE)
if(nrow(mergehash)!=nrow(hash_before)||!all(mergehash$unchanged))stop("Gate4 or Gate4b files changed during correction")

# Source input manifest captures all frozen inputs used by this correction.
inputs<-c(file.path(RAW,"sample_metadata.csv"),file.path(RAW,"series_01_expression_matrix.csv"),file.path(ANN,"GPL96.csv"),file.path(ANN,"GPL97.csv"),
  file.path(PH4,"03_deconvolution_inputs/GSE12021_bulk_expression_for_music.tsv"),file.path(PH4,"03_deconvolution_inputs/GSE12021_metadata_for_deconvolution.csv"),
  file.path(PH4,"01_reference_construction/GSE248453_reference_primary_frozen10_excluded.rds"),
  file.path(BV2,"inputs/reference/reference_B_counts.mtx"),file.path(BV2,"inputs/reference/reference_B_genes.txt"),file.path(BV2,"inputs/reference/reference_B_cell_metadata.csv"),
  file.path(BV2,"inputs/bulk/GSE12021_bulk_expression_for_music.tsv"),
  file.path(BASE,"outputs_reml_hk_robustness/study_level_frozen_inputs.csv"),file.path(BASE,"outputs_reml_hk_robustness/reml_hk_primary_summary.csv"),
  file.path(GATE1,"module_vs_myeloid_dataset_results.csv"),file.path(GATE1,"oa_vs_control_myeloid_results.csv"),file.path(GATE1,"oa_attenuation_dataset_results.csv"),file.path(GATE1,"marker_proxy_sensitivity_results.csv"),
  file.path(BV2,"gate2_Bv2_module_composition.csv"),file.path(BV2,"gate2_Bv2_oa_composition.csv"),file.path(BV2,"gate2_Bv2_attenuation.csv"))
input_manifest<-data.frame(path=normalizePath(inputs,winslash="/"),size_bytes=file.info(inputs)$size,sha256=vapply(inputs,sha256,character(1)),stringsAsFactors=FALSE)
write.csv(input_manifest,file.path(OUT,"source_data/correction_input_manifest.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(OUT,"logs/sessionInfo.txt"))

cat("\nGSE12021 canonical platform: GPL96\n")
cat("unique OA participants:",sum(map$phenotype=="OA"),"\nunique control participants:",sum(map$phenotype=="control"),"\nmodule genes on GPL96:",sum(covA$present),"/10\n")
cat("corrected GSE12021 g:",g120,"95% CI",ci120[1],"to",ci120[2],"\n")
cat("corrected eight-dataset pooled g:",primary_sum$pooled_effect,"95% CI",primary_sum$ci_low,"to",primary_sum$ci_high,"95% PI",primary_sum$pi_low,"to",primary_sum$pi_high,"I2",primary_sum$I2,"\n")
cat("seven-dataset omit-GSE12021 pooled g:",without12021$pooled_effect,"\n")
cat("corrected Reference A pooled rho:",gate1_pool$pooled_effect[gate1_pool$metric=="Reference_A_module_composition_rho"],"\n")
cat("corrected B-v2 pooled rho:",gate2_pool$pooled_effect[gate2_pool$metric=="Bv2_module_composition_rho"],"\n")
cat("corrected Reference A attenuation count:",newA_att_n,"/8\ncorrected B-v2 attenuation count:",newB_att_n,"/8\n")
cat("corrected primary total analysis units:",new_meta_total,"\nscientific impact level:",overall,"\nGate 4 affected: NO\nGate 4b affected: NO\nsafe to render final figures: YES (subject to ledger synchronization and author review)\n")
sink(type="message");sink();close(LOG)
