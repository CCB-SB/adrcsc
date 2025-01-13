library(data.table)
library(Seurat)
library(readxl)
library(lubridate)
library(pbapply)

source("scripts/helper.R")

# cell type (+ subcluster) composition per biogroup with clinical measure, spearman
# 
updrs = as.data.table(read_excel("data/metadata/adrc_updrs_raw.xlsx"))
adrc_nacc = fread("data/metadata/adrc_nacc_data.csv")
adrc_nacc[, visit:=as.numeric(stringr::str_match(redcap_event_name, "visit_(\\d)_arm")[,2])]

adrc_moca = adrc_nacc[, list(moca=as.numeric(gsub("NA", "", paste(c2_mocatots, collapse = '')))), by=c("record_id", "redcap_event_name", "visit")]
adrc_moca[moca == 88, moca:=NA]

adrc_normalcog = adrc_nacc[, list(normalcog=as.numeric(gsub("NA", "", paste(d1_normcog, collapse = '')))), by=c("record_id", "redcap_event_name", "visit")]
adrc_demented = adrc_nacc[, list(demented=as.numeric(gsub("NA", "", paste(d1_demented, collapse = '')))), by=c("record_id", "redcap_event_name", "visit")]

updrs3_cols_start = grep("^updrs_3_1$", colnames(adrc_nacc))
updrs3_cols_end = grep("^updrs_3_18$", colnames(adrc_nacc))

compute_updrs3 = function(x) {
  x[x == 99] = NA
  x[x == 96] = NA
  apply(x, 1, function(e){
    ifelse(mean(is.na(e)) >= 0.3, NA, sum(e, na.rm=T))
  })
}

adrc_nacc[, updrs3:=compute_updrs3(.SD), .SDcols = colnames(adrc_nacc)[updrs3_cols_start:updrs3_cols_end]]

adrc_nacc_for_updrs3 = adrc_nacc[!is.na(updrs3)]
adrc_nacc_for_updrs3 = adrc_nacc_for_updrs3[, date:=mdy(adrc_nacc_for_updrs3$updrs_data_ent_date)]
adrc_nacc_for_updrs3 = adrc_nacc_for_updrs3[order(date)]
adrc_nacc_for_updrs3 = adrc_nacc_for_updrs3[!duplicated(adrc_nacc_for_updrs3, by=c("record_id", "redcap_event_name"), fromLast = TRUE)]

# in case of multiple measures take the one that is nearer to the mean over the visits
adrc_nacc_for_updrs3[, updrs3_mean:=mean(updrs3), by=c("record_id")]
adrc_nacc_for_updrs3 = adrc_nacc_for_updrs3[, updrs3_diff:=abs(updrs3-updrs3_mean)][order(updrs3_diff)]
adrc_nacc_for_updrs3 = adrc_nacc_for_updrs3[!duplicated(adrc_nacc_for_updrs3, by=c("record_id", "redcap_event_name"), fromLast = FALSE)]

adrc_hy = adrc_nacc_for_updrs3[, list(hy=as.numeric(gsub("NA", "", paste(updrs_hoehn_n_yahr, collapse = '')))), by=c("record_id", "redcap_event_name")]

metadata = fread("results/metadata_with_celltype_cluster_proportions_and_new_biomarker_data.csv")
metadata[, `ADRC study ID`:=as.character(`ADRC study ID`)]
metadata = merge(metadata, adrc_nacc_for_updrs3[, c("record_id", "visit", "updrs3")], by.x=c("ADRC study ID", "Visit.x"), by.y=c("record_id", "visit"), all.x=T)
metadata = merge(metadata, adrc_moca[, c("record_id", "visit", "moca")], by.x=c("ADRC study ID", "Visit.x"), by.y=c("record_id", "visit"), all.x=T)

metadata = metadata[!duplicated(Sample)]

# Quanterix values are only present for visit 1
for(quant_c in grep("Quanterix", colnames(metadata), value=T)){
  metadata[Visit.x!=1 & !is.na(metadata[[quant_c]]), sprintf("%s", quant_c):=NA]
}

pbmc = readRDS("results/annotated_celltypes/CompleteObjectAnnotated.rds")

pbmc$biogroup_short = name2short[as.character(pbmc$Diagnosis)]

pbmc$biogroup_short = factor(pbmc@meta.data$biogroup_short, levels=diagnosis_order)

celltype_cluster_stats_per_sample = as.data.table(pbmc@meta.data)[celltype_cluster != "Unknown/doublets"][, list(N=.N), by=c("biogroup_short", "celltype_cluster", "Sample")]
celltype_cluster_stats_per_sample[, prop:=N/sum(N), by="Sample"]

celltype_stats_per_sample = as.data.table(pbmc@meta.data)[celltype != "Unknown/doublets"][, list(N=.N), by=c("biogroup_short", "celltype", "Sample")]
celltype_stats_per_sample[, prop:=N/sum(N), by="Sample"]

for(c in unique(celltype_stats_per_sample$celltype)){
  ct_sub = celltype_stats_per_sample[celltype == c]
  metadata[, sprintf("celltype_%s", c):=ct_sub$prop[match(Sample, ct_sub$Sample)]]
}

for(c in unique(celltype_cluster_stats_per_sample$celltype)){
  ct_sub = celltype_cluster_stats_per_sample[celltype_cluster == c]
  metadata[, sprintf("celltype_cluster_%s", c):=ct_sub$prop[match(Sample, ct_sub$Sample)]]
}

ct_vars = sprintf("celltype_%s", unique(celltype_stats_per_sample$celltype))
ct_clust_vars = sprintf("celltype_cluster_%s", unique(celltype_cluster_stats_per_sample$celltype))

metadata[, `Quanterix3plex.CSF___Ab42/Ab40`:=Quanterix3plex.CSF___Ab42 / Quanterix3plex.CSF___Ab40]
other_vars = c("Age", colnames(metadata)[seq(25,109)], c("updrs3", "moca", "Quanterix3plex.CSF___Ab42/Ab40", "CSF P-Tau181",
                                                         "CSF Total Tau", "CSF AB42", "CSF AB40", "CSF AB42/AB40"))
other_vars = other_vars[!other_vars %in% c("Diagnosis_at_CSF_measurement", "Age_at_CSF_measurement", "ApoE")]


correlations_df = rbindlist(pblapply(c(ct_vars, ct_clust_vars), function(ct){
  res = data.table()
  for(v in other_vars){
    for(d in unique(metadata$Diagnosis)){
      sub = metadata[!is.na(metadata[[ct]]) & !is.na(metadata[[v]]) & Diagnosis == d]
      if(nrow(sub)<=3){
        next;
      }
      pcor = cor.test(sub[[ct]], sub[[v]], method = "pearson")
      scor = cor.test(sub[[ct]], sub[[v]], method = "spearman")
      
      add_covars = c("Age")
      if(length(unique(sub$Sex)) > 1){
        add_covars = c(add_covars, "Sex")
      }
      
      lm_mdl = lm(as.formula(sprintf("`%s` ~ `%s` + %s", v, ct, paste(sprintf('`%s`', add_covars), collapse = ' + '))), sub)
      lm_mdl_p = anova(lm_mdl)[["Pr(>F)"]][1]
      lm_mdl_rsq = summary(lm_mdl)$r.squared
      
      res = rbind(res, data.table(V1=ct, V2=v, Pearson=pcor$estimate, Spearman=scor$estimate, Diagnosis=d, NSamples=nrow(sub),
                                  NPatients=length(unique(sub$SCMD)),
                                  pearson_p=pcor$p.value, spearman_p=scor$p.value,
                                  linear_model_pval=lm_mdl_p,
                                  linear_model_rsq=lm_mdl_rsq,
                                  V1_cv=sd(sub[[ct]])/mean(sub[[ct]]),
                                  V2_cv=sd(sub[[v]])/mean(sub[[v]])))
    }
    sub = metadata[!is.na(metadata[[ct]]) & !is.na(metadata[[v]])]
    if(nrow(sub)<=3){
      next;
    }
    
    pcor = cor.test(sub[[ct]], sub[[v]], method = "pearson")
    scor = cor.test(sub[[ct]], sub[[v]], method = "spearman")
    
    add_covars = c("Age")
    if(length(unique(sub$Sex)) > 1){
      add_covars = c(add_covars, "Sex")
    }
    
    lm_mdl = lm(as.formula(sprintf("`%s` ~ `%s` + %s", v, ct, paste(sprintf('`%s`', add_covars), collapse = ' + '))), sub)
    lm_mdl_p = anova(lm_mdl)[["Pr(>F)"]][1]
    lm_mdl_rsq = summary(lm_mdl)$r.squared
    
    d = "all"
    
    res = rbind(res, data.table(V1=ct, V2=v, Pearson=pcor$estimate, Spearman=scor$estimate, Diagnosis=d, NSamples=nrow(sub),
                                NPatients=length(unique(sub$SCMD)),
                                pearson_p=pcor$p.value, spearman_p=scor$p.value,
                                linear_model_pval=lm_mdl_p,
                                linear_model_rsq=lm_mdl_rsq,
                                V1_cv=sd(sub[[ct]])/mean(sub[[ct]]),
                                V2_cv=sd(sub[[v]])/mean(sub[[v]])))
    
  }
  return(res)
}, cl = 32))

fwrite(correlations_df, "results/correlation/celltype_and_celltype_cluster_correlations_detail.csv")
