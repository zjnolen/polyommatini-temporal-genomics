rule gerp_noN:
    """
    Remove Ns from the GERP scores file. This just makes it faster to load as
    it doesn't have the sites with no info
    """
    input:
        config["gerp_scores"],
    output:
        "results/datasets/{dataset}/analyses/gerp/{dataset}.{ref}_gerps_noNs.tsv.gz",
    conda:
        "../envs/shell.yaml"
    shell:
        """
        zcat {input} | awk '$3 != "N"' | gzip > {output}
        """


rule prep_dtol_gff:
    """
    Sort and index the gene annotation files from DToL
    """
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
    """
    Predict variant effects with VEP using the DToL annotation.
    """
    input:
        bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcf",
        ref="results/ref/{ref}/{ref}.fa",
        gff="results/ref/{ref}/{ref}.gff.gz",
    output:
        out="results/datasets/{dataset}/analyses/vep/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.vep-annotated.txt",
        html="results/datasets/{dataset}/analyses/vep/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.vep-annotated.txt_summary.html",
        warnings="results/datasets/{dataset}/analyses/vep/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.vep-annotated.txt_warnings.txt",
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
    """
    Filter VEP output to include only the picked effect for each variant, i.e.
    the highest impact.
    """
    input:
        vep="results/datasets/{dataset}/analyses/vep/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.vep-annotated.txt",
    output:
        csv="results/datasets/{dataset}/analyses/vep/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.vep-effects.txt",
    conda:
        "../envs/shell.yaml"
    shell:
        r"""
        grep PICK=1 {input.vep} | cut -f2,3,7,14 | sed 's/:/\t/' > {output.csv}
        """


rule bcf2csv:
    """
    Convert BCF to table with only the needed fields, i.e. chrom, pos, ref, alt,
    depth, and genotypes. Will be merged with VEP output.
    """
    input:
        "results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcf",
    output:
        csv="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.csv",
        samples="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.samples",
    conda:
        "../envs/snpsift.yaml"
    shell:
        """
        bcftools view -Ov {input} |
        SnpSift extractFields -s "," -e "." - \
            CHROM POS REF ALT DP "GEN[*].GT" > {output.csv}
        bcftools query -l {input} > {output.samples}
        """


rule calc_mutation_burden:
    """
    Estimates counts of alleles per sample for different variant classes (VEP
    and GERP based). Includes counts when assuming all alternates are
    deleterious or only derived alternates.
    """
    input:
        vep="results/datasets/{dataset}/analyses/vep/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.vep-effects.txt",
        vars="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.csv",
        samples="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.samples",
        pops="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
        anc="results/datasets/{dataset}/analyses/gerp/{dataset}.{ref}_gerps_noNs.tsv.gz",
    output:
        varcounts="results/datasets/{dataset}/analyses/burden/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.varimpacts.tsv",
        varcounts_anc="results/datasets/{dataset}/analyses/burden/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.varimpacts_anc.tsv",
    threads: 24
    conda:
        "../envs/r.yaml"
    resources:
        runtime="6h",
    params:
        gerp_thresh=config["lower_gerp_thresh"],
    script:
        "../scripts/mutation_burden.R"
