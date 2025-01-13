rule diff_exp_muscat:
    input:
        metadata=join(config["results_dir"], "metadata.filtered.csv"),
        obj=join(config["results_dir"], "annotated_celltypes/CompleteObjectAnnotated.rds")
    output:
        sce=join(config["results_dir"], "processed/sce_prep_muscat_first_visit_per_patient_{celltype}.rds"),
        pb=join(config["results_dir"], "processed/pb_sum_prep_muscat_first_visit_per_patient_{celltype}.rds"),
        genes_to_keep=join(config["results_dir"], "processed/genes_to_keep_per_{celltype}.rds"),
        deg_res=expand(join(config["results_dir"], "processed", "ds_pb_{tool}_{{celltype}}.rds"), tool=TOOLS),
        deg_all_col=expand(join(config["results_dir"], "processed", "ds_pb_{tool}_all_col_{{celltype}}.rds"), tool=TOOLS),
    params:
        celltype="{celltype}",
        prefix=join(config["results_dir"], "processed")
    wildcard_constraints: celltype="celltype_cluster" 
    conda: "../envs/single_cell_R_diffexp.yml"
    threads: 128
    script: "../scripts/diff_exp_muscat_celltypes.R"

rule analyse_volcano_dge:
    input:
        metadata=join(config["results_dir"], "metadata.filtered.csv"),
        obj=join(config["results_dir"], "annotated_celltypes", "CompleteObjectAnnotated.rds"),
        deg_all_col=join(config["results_dir"], "processed", "ds_pb_{tool}_all_col_{celltype}.rds"),
        geneid2symbol="data/annotations/geneid2symbol.final.csv"
    output:
        deg_all_col_w_add_info=join(config["results_dir"], "processed", "ds_pb_{tool}_all_col_{celltype}_with_additional_stats.rds"),
        deg_all_col_w_add_info_csv=join(config["results_dir"], "de", "volcano", "ds_pb_{tool}_all_col_{celltype}_with_additional_stats.csv"),
        deg_all_col_w_add_info_sig_p_csv=join(config["results_dir"], "de", "volcano", "ds_pb_{tool}_all_col_{celltype}_with_additional_stats.p_filtered.csv"),
        unfiltered_volcano=join(config["results_dir"], "de", "volcano", "ds_pb_{tool}_volcano_wo_pseudogene_and_unnamed_labels.all.{celltype}.pdf"),
        unfiltered_volcano_unlabeled=join(config["results_dir"], "de", "volcano", "ds_pb_{tool}_volcano_wo_labels.all.{celltype}.pdf"),
        filtered_bio025_c005=join(config["results_dir"], "de", "volcano", "ds_pb_{tool}_volcano_wo_pseudogene_and_unnamed_label.filtered_bio0.25_c0.05.{celltype}.pdf"),
        filtered_bio025_c005_unlabeled=join(config["results_dir"], "de", "volcano", "ds_pb_{tool}_volcano_wo_labels.filtered_bio0.25_c0.05.{celltype}.pdf"),
        filtered_bio05_c005=join(config["results_dir"], "de", "volcano", "ds_pb_{tool}_volcano_wo_pseudogene_and_unnamed_label.filtered_bio0.5_c0.05.{celltype}.pdf"),
        filtered_bio05_c005_unlabeled=join(config["results_dir"], "de", "volcano", "ds_pb_{tool}_volcano_wo_labels.filtered_bio0.5_c0.05.{celltype}.pdf"),
    params:
        celltype_annot="{celltype}"
    #wildcard_constraints: celltype= "|".join(celltype_annotations.keys()) + "|" + "|".join("{}_w_batch".format(c) for c in celltype_annotations.keys()) + "|" +  "|".join("{}_w_batch_and_sex".format(c) for c in celltype_annotations.keys())
    wildcard_constraints: celltype="celltype_cluster|celltype"
    conda: "../envs/single_cell_R_diffexp.yml"
    threads: 32
    script: "../scripts/analyse_volcano_dge.R"


rule analyse_general_dge:
    input:
        deg_all_col_w_add_info=join(config["results_dir"], "processed", "ds_pb_{tool}_all_col_{celltype}_with_additional_stats.rds"),
        obj=join(config["results_dir"], "annotated_celltypes", "CompleteObjectAnnotated.rds"),
    output:
        all_comp_heat=join(config["results_dir"], "de", "heatmap", "ds_pb_{tool}_all_comparisons_{celltype}.all.pdf"),
        all_comp_heat_bio025_c005=join(config["results_dir"], "de", "heatmap", "ds_pb_{tool}_all_comparisons_{celltype}.filtered_bio0.25_c0.05.pdf"),
        all_comp_heat_bio05_c005=join(config["results_dir"], "de", "heatmap", "ds_pb_{tool}_all_comparisons_{celltype}.filtered_bio0.5_c0.05.pdf"),
        up_comp_heat=join(config["results_dir"], "de", "heatmap", "ds_pb_{tool}_up_comparisons_{celltype}.all.pdf"),
        up_comp_heat_bio025_c005=join(config["results_dir"], "de", "heatmap", "ds_pb_{tool}_up_comparisons_{celltype}.filtered_bio0.25_c0.05.pdf"),
        up_comp_heat_bio05_c005=join(config["results_dir"], "de", "heatmap", "ds_pb_{tool}_up_comparisons_{celltype}.filtered_bio0.5_c0.05.pdf"),
        down_comp_heat=join(config["results_dir"], "de", "heatmap", "ds_pb_{tool}_down_comparisons_{celltype}.all.pdf"),
        down_comp_heat_bio025_c005=join(config["results_dir"], "de", "heatmap", "ds_pb_{tool}_down_comparisons_{celltype}.filtered_bio0.25_c0.05.pdf"),
        down_comp_heat_bio05_c005=join(config["results_dir"], "de", "heatmap", "ds_pb_{tool}_down_comparisons_{celltype}.filtered_bio0.5_c0.05.pdf"),
        all_comp_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_all_comparisons_{celltype}.all.pdf"),
        all_comp_upset_bio025_c005=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_all_comparisons_{celltype}.filtered_bio0.25_c0.05.pdf"),
        all_comp_upset_bio05_c005=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_all_comparisons_{celltype}.filtered_bio0.5_c0.05.pdf"),
        up_comp_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_up_comparisons_{celltype}.all.pdf"),
        up_comp_upset_bio025_c005=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_up_comparisons_{celltype}.filtered_bio0.25_c0.05.pdf"),
        up_comp_upset_bio05_c005=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_up_comparisons_{celltype}.filtered_bio0.5_c0.05.pdf"),
        down_comp_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_down_comparisons_{celltype}.all.pdf"),
        down_comp_upset_bio025_c005=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_down_comparisons_{celltype}.filtered_bio0.25_c0.05.pdf"),
        down_comp_upset_bio05_c005=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_down_comparisons_{celltype}.filtered_bio0.5_c0.05.pdf"),
        disease_comp_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_all_comparisons_{celltype}.all.diagnosis_only.pdf"),
        disease_comp_bio025_c005_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_all_comparisons_{celltype}.filtered_bio0.25_c0.05.diagnosis_only.pdf"),
        disease_comp_bio05_c005_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_all_comparisons_{celltype}.filtered_bio0.5_c0.05.diagnosis_only.pdf"),
        disease_up_comp_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_up_comparisons_{celltype}.all.diagnosis_only.pdf"),
        disease_up_comp_bio025_c005_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_up_comparisons_{celltype}.filtered_bio0.25_c0.05.diagnosis_only.pdf"),
        disease_up_comp_bio05_c005_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_up_comparisons_{celltype}.filtered_bio0.5_c0.05.diagnosis_only.pdf"),
        disease_down_comp_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_down_comparisons_{celltype}.all.diagnosis_only.pdf"),
        disease_down_comp_bio025_c005_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_down_comparisons_{celltype}.filtered_bio0.25_c0.05.diagnosis_only.pdf"),
        disease_down_comp_bio05_c005_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_down_comparisons_{celltype}.filtered_bio0.5_c0.05.diagnosis_only.pdf"),
        celltype_comp_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_all_comparisons_{celltype}.all.celltype_only.pdf"),
        celltype_comp_bio025_c005_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_all_comparisons_{celltype}.filtered_bio0.25_c0.05.celltype_only.pdf"),
        celltype_comp_bio05_c005_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_all_comparisons_{celltype}.filtered_bio0.5_c0.05.celltype_only.pdf"),
        celltype_up_comp_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_up_comparisons_{celltype}.all.celltype_only.pdf"),
        celltype_up_comp_bio025_c005_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_up_comparisons_{celltype}.filtered_bio0.25_c0.05.celltype_only.pdf"),
        celltype_up_comp_bio05_c005_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_up_comparisons_{celltype}.filtered_bio0.5_c0.05.celltype_only.pdf"),
        celltype_down_comp_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_down_comparisons_{celltype}.all.celltype_only.pdf"),
        celltype_down_comp_bio025_c005_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_down_comparisons_{celltype}.filtered_bio0.25_c0.05.celltype_only.pdf"),
        celltype_down_comp_bio05_c005_upset=join(config["results_dir"], "de", "upset", "ds_pb_{tool}_down_comparisons_{celltype}.filtered_bio0.5_c0.05.celltype_only.pdf"),
        deg_vs_cells=join(config["results_dir"], "de", "correlation", "ds_pb_{tool}_all_comparisons_deg_vs_cells_{celltype}.all.pdf"),
        deg_vs_cells_bio025_c005=join(config["results_dir"], "de", "correlation", "ds_pb_{tool}_all_comparisons_deg_vs_cells_{celltype}.filtered_bio0.25_c0.05.pdf"),
        deg_vs_cells_bio05_c005=join(config["results_dir"], "de", "correlation", "ds_pb_{tool}_all_comparisons_deg_vs_cells_{celltype}.filtered_bio0.5_c0.05.pdf"),
        deregulation_dir=join(config["results_dir"], "de", "deregulation_dir", "ds_pb_{tool}_all_comparisons_{celltype}.total.all.pdf"),
        deregulation_dir_bio025_c005=join(config["results_dir"], "de", "deregulation_dir", "ds_pb_{tool}_all_comparisons_{celltype}.total.filtered_bio0.25_c0.05.pdf"),
        deregulation_dir_bio05_c005=join(config["results_dir"], "de", "deregulation_dir", "ds_pb_{tool}_all_comparisons_{celltype}.total.filtered_bio0.5_c0.05.pdf"),
        deregulation_dir_perc=join(config["results_dir"], "de", "deregulation_dir", "ds_pb_{tool}_all_comparisons_{celltype}.percent.all.pdf"),
        deregulation_dir_bio025_c005_perc=join(config["results_dir"], "de", "deregulation_dir", "ds_pb_{tool}_all_comparisons_{celltype}.percent.filtered_bio0.25_c0.05.pdf"),
        deregulation_dir_bio05_c005_perc=join(config["results_dir"], "de", "deregulation_dir", "ds_pb_{tool}_all_comparisons_{celltype}.percent.filtered_bio0.5_c0.05.pdf")
    params:
        celltype_annot="{celltype}"
    wildcard_constraints: celltype="celltype_cluster|celltype"
    conda: "../envs/rstudio_env.yml"
    threads: 1
    script: "../scripts/plot_deg.R"
