rule bcftools_mpileup_call_individual:
    """
    Call genotypes per sample in filtered sites regions.
    """
    input:
        alignments="results/datasets/{dataset}/bams/{population}.{ref}{dp}.bam",
        bai="results/datasets/{dataset}/bams/{population}.{ref}{dp}.bam.bai",
        ref="results/ref/{ref}/{ref}.fa",
        index="results/ref/{ref}/{ref}.fa.fai",
        regions="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.bed",
    output:
        bcf=temp(
            "results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.bcf"
        ),
        idx=temp(
            "results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.bcf.csi"
        ),
        stats="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.bcf.stats",
    conda:
        "../envs/bcftools.yaml"
    params:
        baseq=config["baseQ"],
        mapq=config["mapQ"],
    log:
        "logs/{dataset}/bcftools/mpileup/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    threads: 2
    resources:
        runtime=720,
    shell:
        """
        bcftools mpileup --threads {threads} -f {input.ref} -R {input.regions} -Ou \
            -B --min-MQ {params.mapq} --min-BQ {params.baseq} \
            -a "FORMAT/AD,FORMAT/DP,INFO/AD" {input.alignments} |
        bcftools call -m -f GQ,GP -o {output.bcf}
        bcftools index -o {output.idx} {output.bcf}
        bcftools stats {output.bcf} > {output.stats}
        """


rule bcftools_filter:
    """
    Filter out indels, SNPs within 5 bp of indels, and set genotypes with a
    depth less than a specified depth to missing. Heterozygotes are set to
    missing if they don't have an allelic balance between 0.2 and 0.8.
    """
    input:
        bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.bcf",
        idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.bcf.csi",
    output:
        bcf=temp(
            "results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}.bcf"
        ),
        idx=temp(
            "results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}.bcf.csi"
        ),
        stats="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}.bcf.stats",
    conda:
        "../envs/bcftools.yaml"
    threads: 2
    resources:
        runtime=720,
    shell:
        """
        bcftools filter -g 5 -Ou {input.bcf} | \
            bcftools view -V indels -Ou | \
            bcftools +setGT -Ou -- -t q -n . -i"FMT/DP<{wildcards.mindp}" | \
            bcftools +setGT -Ou -- -t q -n . -i'GT="het" & (FMT/AD[:0]/FMT/DP < 0.2 | FMT/AD[:0]/FMT/DP > 0.8 | FMT/AD[:1]/FMT/DP < 0.2 | FMT/AD[:1]/FMT/DP > 0.8)' | \
            bcftools +fill-tags -Ob -- -t all > {output.bcf}
        bcftools index -o {output.idx} {output.bcf}
        bcftools stats {output.bcf} > {output.stats}
        """


rule bcftools_filter_nomiss:
    """
    Remove the now missing positions from the per sample bcf.
    """
    input:
        bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}.bcf",
        idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}.bcf.csi",
    output:
        bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}.nomiss.bcf",
        idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}.nomiss.bcf.csi",
        stats="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}.nomiss.bcf.stats",
    conda:
        "../envs/bcftools.yaml"
    threads: 2
    resources:
        runtime=720,
    shell:
        """
        bcftools view -i'GT!="./."' -Ob {input.bcf} > {output.bcf}
        bcftools index -o {output.idx} {output.bcf}
        bcftools stats {output.bcf} > {output.stats}
        """


def bcfs(wildcards):
    dp = wildcards.dp
    samples = [
        sample
        for sample in angsd.samples.index.tolist()
        if sample not in config["calling_drop"]
    ]
    return expand(
        "results/datasets/{{dataset}}/bcfs/{{dataset}}.{{ref}}_{population}{{dp}}_{{sites}}-filts.filtered_mindp{{mindp}}.nomiss.bcf",
        population=samples,
    )


def bcf_stats(wildcards):
    dp = wildcards.dp
    samples = [
        sample
        for sample in angsd.samples.index.tolist()
        if sample not in config["calling_drop"]
    ]
    return expand(
        "results/datasets/{{dataset}}/bcfs/{{dataset}}.{{ref}}_{population}{{dp}}_{{sites}}-filts.filtered_mindp{{mindp}}.nomiss.bcf.stats",
        population=samples,
    )


rule multiqc_bcf_filt_nomiss:
    """
    Generate multiqc with filtered BCF stats for each individual.
    """
    input:
        bcf_stats,
    output:
        "results/datasets/{dataset}/qc/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}.nomiss.bcf.multiqc.html",
    params:
        extra="--verbose",
        use_input_files_only=True,
    wrapper:
        "v3.7.0/bio/multiqc"


rule bcftools_merge_all:
    """
    Generate dataset-wide BCF from per sample BCFs. Sites are filtered to
    only biallelic sites with at least one minor allele present in the
    dataset. Sites are also filtered by a missing data threshold.
    """
    input:
        bcfs,
    output:
        bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.trans.fmiss{miss}.bcf",
        idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.trans.fmiss{miss}.bcf.csi",
        stats="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.trans.fmiss{miss}.bcf.stats",
    conda:
        "../envs/bcftools.yaml"
    threads: 2
    resources:
        runtime=720,
    shell:
        """
        bcftools merge -Ou {input} |
            bcftools +fill-tags -Ou -- -t all |
            bcftools view -m2 -M2 -v snps -i 'MAF > 0' -Ou |
            bcftools view -i 'F_MISSING <= {wildcards.miss}' -Ob > {output.bcf}
        bcftools index -o {output.idx} {output.bcf}
        bcftools stats {output.bcf} > {output.stats}
        # bcftools view -m2 -M2 -v snps -e 'MAF < 0.05' -Ou |
        """


rule bcftools_remove_trans:
    input:
        bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.trans.fmiss{miss}.bcf",
        idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.trans.fmiss{miss}.bcf.csi",
    output:
        bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.notrans.fmiss{miss}.bcf",
        idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.notrans.fmiss{miss}.bcf.csi",
        stats="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.notrans.fmiss{miss}.bcf.stats",
    conda:
        "../envs/bcftools.yaml"
    threads: 2
    resources:
        runtime=720,
    shell:
        r"""
        bcftools view -Ov {input.bcf} |
            awk -F '\t' '!(($4 == "A" && $5 == "G") || ($4 == "G" && $5 == "A") || ($4 == "C" && $5 == "T") || ($4 == "T" && $5 == "C"))' |
            bcftools view -Ob > {output.bcf}
        bcftools index -o {output.idx} {output.bcf}
        bcftools stats {output.bcf} > {output.stats}
        """
