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


# localrules: bcftools_list
#
# rule bcftools_list:
#     """
#     Create BED file of chromosome groups to partition variant calling in
#     BCFtools with.
#     """
#     input:
#         bed="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.bed",
#         reg="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
#     output:
#         bed="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_chunk{chunk}_{sites}-filts.bed",
#     conda:
#         "../envs/shell.yaml"
#     shell:
#         """
#         > {output.bed}
#
#         while read chr; do
#             awk -v chr=$chr '$1 == chr' {input.bed} >> {output.bed}
#         done < {input.reg}
#         """
#
#
# rule bcftools_mpileup_call:
#     """
#     Call genotypes across the dataset in filtered sites regions.
#     """
#     input:
#         alignments=angsd.get_bamlist_bams,
#         bai=angsd.get_bamlist_bais,
#         ref="results/ref/{ref}/{ref}.fa",
#         index="results/ref/{ref}/{ref}.fa.fai",
#         regions="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_chunk{chunk}_{sites}-filts.bed",
#     output:
#         bcf=temp(
#             "results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.bcf"
#         ),
#         idx=temp(
#             "results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.bcf.csi"
#         ),
#     wildcard_constraints:
#         dp=".{0}",
#     conda:
#         "../envs/bcftools.yaml"
#     params:
#         baseq=config["baseQ"],
#         mapq=config["mapQ"],
#     log:
#         "logs/{dataset}/bcftools/mpileup/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log",
#     threads: 2
#     resources:
#         runtime=1440,
#     shell:
#         """
#         bcftools mpileup --threads {threads} -f {input.ref} -R {input.regions} -Ou \
#             -B --min-MQ {params.mapq} --min-BQ {params.baseq} \
#             -a "FORMAT/AD,FORMAT/DP,INFO/AD" {input.alignments} |
#         bcftools call -m -f GQ,GP -o {output.bcf}
#
#         bcftools index -o {output.idx} {output.bcf}
#         """
#
#
# rule bcf_concat:
#     input:
#         expand(
#             "results/datasets/{{dataset}}/bcfs/{{dataset}}.{{ref}}_{{population}}_chunk{chunk}_{{sites}}-filts.bcf",
#             chunk=angsd.chunklist,
#         ),
#     output:
#         bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.bcf",
#         idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.bcf.csi",
#         stats="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.bcf.stats",
#     conda:
#         "../envs/bcftools.yaml"
#     threads: 2
#     resources:
#         runtime=720,
#     shell:
#         """
#         bcftools concat -Ob {input} > {output.bcf}
#         bcftools index -o {output.idx} {output.bcf}
#         bcftools stats {output.bcf} > {output.stats}
#         """
#
#
# rule bcftools_filter:
#     """
#     Filter out indels, SNPs within 5 bp of indels, and set genotypes with a
#     depth less than 4 or genotype quality less than 30 to missing. In
#     testing, GQ>=30 is as stringent a filter as an allelic balance between
#     0.2 and 0.8 to call a heterozygote, so we don't filter for allelic balance.
#     """
#     input:
#         bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.bcf",
#         idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.bcf.csi",
#     output:
#         bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.filtered-biallelic.trans.fmiss{miss}.bcf",
#         idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.filtered-biallelic.trans.fmiss{miss}.bcf.csi",
#         stats="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.filtered-biallelic.trans.fmiss{miss}.bcf.stats",
#     conda:
#         "../envs/bcftools.yaml"
#     threads: 2
#     resources:
#         runtime=720,
#     shell:
#         """
#         bcftools filter -g 5 -Ou {input.bcf} | \
#             bcftools view -V indels -Ou | \
#             bcftools +setGT -Ou -- -t q -n . -i"FMT/DP<4 | FMT/GQ<30" | \
#             bcftools +fill-tags -Ou -- -t all |
#             bcftools view -m2 -M2 -v snps -e 'MAF < 0.05' -Ou |
#             bcftools view -i 'F_MISSING <= {wildcards.miss}' -Ob > {output.bcf}
#         bcftools index -o {output.idx} {output.bcf}
#         bcftools stats {output.bcf} > {output.stats}
#         """


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


# rule multiqc_bcf_unfilt:
#     """
#     Generate multiqc with unfiltered BCF stats for each individual.
#     """
#     input:
#         expand(
#             "results/datasets/{{dataset}}/bcfs/{{dataset}}.{{ref}}_{population}_{{sites}}-filts.bcf.stats",
#             population=angsd.samples.index,
#         ),
#     output:
#         "results/datasets/{dataset}/qc/bcfs/{dataset}.{ref}_all_{sites}-filts.bcf.multiqc.html",
#     params:
#         extra="--verbose",
#         use_input_files_only=True,
#     wrapper:
#         "v3.7.0/bio/multiqc"





# rule bcftools_merge_pops:
#     """
#     Generate per population BCF from per sample BCFs. Sites are filtered to
#     only biallelic sites with at least one minor allele present in the
#     population. Sites are also filtered by a missing data threshold.
#     """
#     input:
#         bcf=lambda w: expand(
#             "results/datasets/{{dataset}}/bcfs/{{dataset}}.{{ref}}_{population}_{{sites}}-filts.filtered.nomiss.bcf",
#             population=angsd.samples.index[
#                 angsd.samples["population"] == w.population
#             ],
#         ),
#         idx=lambda w: expand(
#             "results/datasets/{{dataset}}/bcfs/{{dataset}}.{{ref}}_{population}_{{sites}}-filts.filtered.nomiss.bcf.csi",
#             population=angsd.samples.index[
#                 angsd.samples["population"] == w.population
#             ],
#         ),
#     output:
#         bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.filtered-biallelic.trans.fmiss{miss}.bcf",
#         idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.filtered-biallelic.trans.fmiss{miss}.bcf.csi",
#         stats="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.filtered-biallelic.trans.fmiss{miss}.bcf.stats",
#     wildcard_constraints:
#         population="|".join(angsd.pop_list),
#     conda:
#         "../envs/bcftools.yaml"
#     threads: 2
#     resources:
#         runtime=720,
#     shell:
#         """
#         bcftools merge -Ou {input.bcf} |
#             bcftools +fill-tags -Ou -- -t all |
#             bcftools view -m2 -M2 -v snps -e 'MAC == 0' -Ou |
#             bcftools view -i 'F_MISSING <= {wildcards.miss}' -Ob > {output.bcf}
#         bcftools index -o {output.idx} {output.bcf}
#         bcftools stats {output.bcf} > {output.stats}
#         """

def bcfs(wildcards):
    dp = wildcards.dp
    samples=[sample for sample in angsd.samples.index.tolist() if sample not in config["calling_drop"]]
    return expand(
        "results/datasets/{{dataset}}/bcfs/{{dataset}}.{{ref}}_{population}{{dp}}_{{sites}}-filts.filtered_mindp{{mindp}}.nomiss.bcf",
        population=samples
    )


def bcf_stats(wildcards):
    dp = wildcards.dp
    samples=[sample for sample in angsd.samples.index.tolist() if sample not in config["calling_drop"]]
    return expand(
        "results/datasets/{{dataset}}/bcfs/{{dataset}}.{{ref}}_{population}{{dp}}_{{sites}}-filts.filtered_mindp{{mindp}}.nomiss.bcf.stats",
        population=samples
    )


rule multiqc_bcf_filt_nomiss:
    """
    Generate multiqc with filtered BCF stats for each individual.
    """
    input:
        bcf_stats
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
        bcfs
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


# rule multiqc_bcf_populations:
#     """
#     Generate multiqc with all population BCF stats.
#     """
#     input:
#         expand(
#             "results/datasets/{{dataset}}/bcfs/{{dataset}}.{{ref}}_{population}_{{sites}}-filts.filtered-biallelic.{{trans}}.fmiss{{miss}}.bcf.stats",
#             population=angsd.pop_list,
#         ),
#     output:
#         "results/datasets/{dataset}/qc/bcfs/{dataset}.{ref}_populations_{sites}-filts.filtered-biallelic.{trans}.fmiss{miss}.bcf.multiqc.html",
#     params:
#         extra="--verbose",
#         use_input_files_only=True,
#     wrapper:
#         "v3.7.0/bio/multiqc"


# rule bcftools_merge_groups:
#     input:
#         bcf=expand(
#             "results/datasets/{{dataset}}/bcfs/{{dataset}}.{{ref}}_{group}_{{sites}}-filts.bcf",
#             group=["modern", "historical"],
#         ),
#     output:
#         bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all_{sites}-filts_fmiss{miss}.bcf",
#         idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all_{sites}-filts_fmiss{miss}.bcf.csi",
#         stats="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all_{sites}-filts_fmiss{miss}.bcf.stats",
#     conda:
#         "../envs/bcftools.yaml"
#     threads: 2
#     resources:
#         runtime=720,
#     shell:
#         """
#         bcftools merge -Ou {input.bcf} |
#             bcftools +fill-tags -Ou -- -t all |
#             bcftools view -i 'F_MISSING <= {wildcards.miss}' -Ou > {output.bcf}
#         bcftools index -o {output.idx} {output.bcf}
#         bcftools stats {output.bcf} > {output.stats}
#         """


# rule bcftools_bcf2bed:
#     input:
#         bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all_{sites}-filts_fmiss{miss}.bcf",
#     output:
#         bed="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all_{sites}-filts_fmiss{miss}.bed",
#     conda:
#         "../envs/bcftools.yaml"
#     threads: 2
#     resources:
#         runtime=720,
#     shell:
#         """
#         bcftools view {input.bcf} | grep -v "^#" | awk '{{print $1"\t"$2-1"\t"$2}}' > {output.bed}
#         """
# rule bcftools_subset_sample_bcf:
#     input:
#         bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.filtered.nomiss.bcf",
#         idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.filtered.nomiss.bcf.csi",
#         bed="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all_{sites}-filts_fmiss{miss}.bed",
#     output:
#         bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.filtered.fmiss{miss}.bcf",
#         idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.filtered.fmiss{miss}.bcf.csi",
#         stats="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}_{sites}-filts.filtered.fmiss{miss}.bcf.stats",
#     conda:
#         "../envs/bcftools.yaml"
#     threads: 2
#     resources:
#         runtime=720,
#     shell:
#         """
#         bcftools view -Ob -T {input.bed} {input.bcf} > {output.bcf}
#         bcftools index -o {output.idx} {output.bcf}
#         bcftools stats {output.bcf} > {output.stats}
#         """
# rule multiqc_bcf_filt_fmiss:
#     input:
#         expand(
#             "results/datasets/{{dataset}}/bcfs/{{dataset}}.{{ref}}_{population}_{{sites}}-filts.filtered.fmiss{{miss}}.bcf.stats",
#             population=angsd.samples.index,
#         ),
#     output:
#         "results/datasets/{dataset}/qc/bcfs/{dataset}.{ref}_all_{sites}-filts.fmiss{miss}.bcf.multiqc.html",
#     params:
#         extra="--verbose",
#         use_input_files_only=True,
#     wrapper:
#         "v3.7.0/bio/multiqc"


rule bcftools_roh_indiv_bcf:
    """
    Call genotypes per sample in filtered sites regions.
    """
    input:
        alignments="results/datasets/{dataset}/bams/{sample}.{ref}{dp}.bam",
        bai="results/datasets/{dataset}/bams/{sample}.{ref}{dp}.bam.bai",
        ref="results/ref/{ref}/{ref}.fa",
        index="results/ref/{ref}/{ref}.fa.fai",
        regions="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.bed",
    output:
        bcf="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{sample}{dp}_{sites}-filts.filtered.notrans.bcf",
        roh="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{sample}{dp}_{sites}-filts.filtered.notrans.roh",
        sites="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{sample}{dp}_{sites}-filts.filtered.notrans.sites.roh",
        regs="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{sample}{dp}_{sites}-filts.filtered.notrans.regs.roh",
    conda:
        "../envs/bcftools.yaml"
    params:
        baseq=config["baseQ"],
        mapq=config["mapQ"],
        recrate=config["recrate"],
    log:
        "logs/{dataset}/bcftools/roh/{dataset}.{ref}_{sample}{dp}_{sites}-filts.log",
    threads: 2
    resources:
        runtime="2d",
    shell:
        """
        bcftools mpileup --threads {threads} -f {input.ref} -R {input.regions} -Ou \
            -B --min-MQ {params.mapq} --min-BQ {params.baseq} \
            -a "FORMAT/AD,FORMAT/DP,INFO/AD" {input.alignments} |
            bcftools call -mv -f GQ,GP -Ou |
            bcftools filter -g 5 -Ou | \
            bcftools view -V indels -Ou | \
            bcftools +setGT -Ou -- -t q -n . -i"FMT/DP<3 | FMT/GQ<30" | \
            bcftools +fill-tags -Ov -- -t all |
            awk -F '\t' '!(($4 == "A" && $5 == "G") || ($4 == "G" && $5 == "A") || ($4 == "C" && $5 == "T") || ($4 == "T" && $5 == "C"))' |
            bcftools view -Ob > {output.bcf}
        bcftools roh -G30 -M {params.recrate} --AF-dflt 0.4 -o {output.roh} {output.bcf}
        awk '$1=="ST"' {output.roh} > {output.sites}
        awk '$1=="RG"' {output.roh} > {output.regs}
        """

def rohs(wildcards):
    dp = wildcards.dp
    samples=[sample for sample in angsd.samples.index.tolist() if sample not in config["calling_drop"]]
    return expand(
        "results/datasets/{{dataset}}/analyses/roh/bcftools/{{dataset}}.{{ref}}_{sample}{{dp}}_{{sites}}-filts.filtered.notrans.regs.roh",
        sample=samples
    )


rule merge_bcftools_roh_indiv:
    input:
        rohs
    output:
        "results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.filtered.notrans.regs.roh"
    conda:
        "../envs/shell.yaml"
    shell:
        """
        cat {input} > {output}
        """

rule bcftools_roh_indiv_plot:
    input:
        roh="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.filtered.notrans.regs.roh",
        inds="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
        autos=angsd.get_auto_sum
    output:
        barplot=report(
            "results/datasets/{dataset}/plots/inbreeding/{dataset}.{ref}_all{dp}_{sites}-filts.filtered.notrans.bcftools.froh_bins.svg",
            category="06 Inbreeding",
            labels=lambda w: {
                "Filter": "{sites}",
                **angsd.dp_report(w),
                "Type": "Froh Bins Barplot (bcftools, ind bcf)",
            },
        ),
        scatter=report(
            "results/datasets/{dataset}/plots/inbreeding/{dataset}.{ref}_all{dp}_{sites}-filts.filtered.notrans.bcftools.cumroh_nroh.svg",
            category="06 Inbreeding",
            labels=lambda w: {
                "Filter": "{sites}",
                **angsd.dp_report(w),
                "Type": "Nroh ~ Lroh Scatterplot (bcftools, ind bcf)",
            },
        ),
        froh="results/datasets/{dataset}/plots/inbreeding/{dataset}.{ref}_all{dp}_{sites}-filts.filtered.notrans.bcftools.ind_froh.tsv",
    log:
        "logs/{dataset}/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.filtered.notrans.bcftools_plot.log",
    conda:
        "../envs/r.yaml"
    params:
        bins=config["params"]["ngsf-hmm"]["roh_bins"],
        minroh=config["params"]["ngsf-hmm"]["min_roh_length"],
        outpre=lambda w, output: output["barplot"].removesuffix(".froh_bins.svg"),
        roh_phred=config["min_phred_roh"],
    script:
        "../scripts/plot_Froh.R"