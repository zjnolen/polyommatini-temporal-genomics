rule prep_dtol_gff:
    input:
        config["gff"],
    output:
        "results/ref/{ref}/{ref}.gff.gz",
    conda:
        "../envs/vep.yaml"
    shell:
        """
        zcat {input} | sort -k1,1d -k4,4n | bgzip > {output}
        tabix {output}
        """


rule vep_annotate_vars:
    input:
        bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.bcf",
        ref="results/ref/{ref}/{ref}.fa",
        gff="results/ref/{ref}/{ref}.gff.gz",
    output:
        out="results/datasets/{dataset}/analyses/vep/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.vep-annotated.txt",
        html="results/datasets/{dataset}/analyses/vep/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.vep-annotated.txt_summary.html",
        warnings="results/datasets/{dataset}/analyses/vep/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.vep-annotated.txt_warnings.txt",
    conda:
        "../envs/vep.yaml"
    resources:
        runtime="4d",
    shell:
        """
        bcftools view -Ov {input.bcf} | \
        vep --fasta {input.ref} --gff {input.gff} --format vcf \
            --flag_pick --force_overwrite -o {output.out}
        """


rule vep_picked_effect:
    input:
        vep="results/datasets/{dataset}/analyses/vep/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.vep-annotated.txt",
    output:
        csv="results/datasets/{dataset}/analyses/vep/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.vep-effects.txt",
    conda:
        "../envs/shell.yaml"
    shell:
        r"""
        grep PICK=1 {input.vep} | cut -f2,3,7,14 | sed 's/:/\t/' > {output.csv}
        """


rule calc_vep_load:
    """
    Estimates normalized counts of impact classes and rxy for impact classes per
    time period.
    """
    input:
        vep="results/datasets/{dataset}/analyses/vep/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.vep-effects.txt",
        vars="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.csv",
        samples="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.samples",
        pops="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
    output:
        varcounts="results/datasets/{dataset}/analyses/vep/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.varimpacts.tsv",
        varcounts_boots="results/datasets/{dataset}/analyses/vep/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.varimpacts_boots.tsv",
        varcounts_nomiss="results/datasets/{dataset}/analyses/vep/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.varimpacts_nomiss.tsv",
    threads: lambda w, attempt: 2 * attempt
    conda:
        "../envs/r.yaml"
    resources:
        runtime="6h",
    script:
        "../scripts/genetic_load_vep.R"
