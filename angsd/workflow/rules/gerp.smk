rule filter_gerps:
    input:
        config["gerp_scores"],
    output:
        "results/datasets/{dataset}/analyses/gerp/{dataset}.{ref}_gerps.csv",
    conda:
        "../envs/shell.yaml"
    params:
        thresh=config["lower_gerp_thresh"],
    shell:
        """
        zcat {input} | awk '$4 >= {params.thresh} && $3 != "N"' > {output}
        """


rule bcf2csv:
    input:
        "results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.bcf",
    output:
        csv="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.csv",
        samples="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.samples",
    conda:
        "../envs/snpsift.yaml"
    shell:
        """
        bcftools view -Ov {input} |
        SnpSift extractFields -s "," -e "." - \
            CHROM POS REF ALT DP "GEN[*].GT" > {output.csv}
        bcftools query -l {input} > {output.samples}
        """


rule calculate_relative_load:
    input:
        csv="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.csv",
        gerp="results/datasets/{dataset}/analyses/gerp/{dataset}.{ref}_gerps.csv",
        samples="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.samples",
        samplemeta="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
    output:
        load="results/datasets/{dataset}/analyses/gerp/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.relative-load-gerp.csv",
        load_nomiss="results/datasets/{dataset}/analyses/gerp/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.relative-load-gerp.nomiss.csv",
        load_boot="results/datasets/{dataset}/analyses/gerp/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.relative-load-gerp.boots.csv",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/relative_gen_load_gerp.R"
