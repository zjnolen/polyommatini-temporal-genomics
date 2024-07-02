rule bcftools_roh:
    """
    Estimate runs of homozygosity with BCFtools. We call and filter genotypes
    per individual. We call and filter genotypes per sample, then merge all
    samples of the same species, removing transitions and invariant sites. We
    then run bcftools roh, ignoring homozygous reference genotypes and setting
    a default allele frequency of 0.4, which mirrors the way the analysis works
    on single sample vcfs, but with the benefit of removing transitions globally
    and ensuring fixed alternate variants are excluded (rather than considered
    to be segregating in the sample set, which they are likely to not be and is
    important for the species we look at that don't have species specific refs).
    This should make the method comparable enough across species despite having
    different sample sizes and having different substructure that makes
    calculating true allele frequencies challenging.
    """
    input:
        "results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.bcf",
    output:
        roh="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.roh",
        sites="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.sites.roh",
        regs="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.regs.roh",
    log:
        "logs/{dataset}/bcftools/roh/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.log",
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards, attempt: attempt
    params:
        recrate=config["recrate"],
    shell:
        """
        (bcftools roh -G30 --ignore-homref -M {params.recrate} --AF-dflt 0.4 \
            -o {output.roh} {input}
        awk '$1=="ST"' {output.roh} > {output.sites}
        awk '$1=="RG"' {output.roh} > {output.regs}) 2> {log}
        """


rule bcftools_roh_plot:
    input:
        roh="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.regs.roh",
        inds="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
        autos=angsd.get_auto_sum,
    output:
        barplot=report(
            "results/datasets/{dataset}/plots/inbreeding/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.bcftools.froh_bins.svg",
            category="06 Inbreeding",
            labels=lambda w: {
                "Filter": "{sites}",
                **angsd.dp_report(w),
                "minCallDP": "{mindp}",
                "Type": "Froh Bins Barplot (bcftools)",
            },
        ),
        scatter=report(
            "results/datasets/{dataset}/plots/inbreeding/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.bcftools.cumroh_nroh.svg",
            category="06 Inbreeding",
            labels=lambda w: {
                "Filter": "{sites}",
                **angsd.dp_report(w),
                "minCallDP": "{mindp}",
                "Type": "Nroh ~ Lroh Scatterplot (bcftools)",
            },
        ),
        froh="results/datasets/{dataset}/plots/inbreeding/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.bcftools.ind_froh.tsv",
    log:
        "logs/{dataset}/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts_filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}_plot.log",
    conda:
        "../envs/r.yaml"
    params:
        bins=config["params"]["ngsf-hmm"]["roh_bins"],
        minroh=config["params"]["ngsf-hmm"]["min_roh_length"],
        outpre=lambda w, output: output["barplot"].removesuffix(".froh_bins.svg"),
        roh_phred=config["min_phred_roh"],
    script:
        "../scripts/plot_Froh.R"


# localrules:
#     combine_roh,
#     gatk_list,
#     bcftools_list,
#     bcf2plink,

# import math

# rule angsd_doBcf_likes:
#     """
#     Generate per population BCF files from ANGSD containing the same genotype
#     likelihoods as used in the other GL analyses.
#     """
#     input:
#         bam="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
#         bams=angsd.get_bamlist_bams,
#         bais=angsd.get_bamlist_bais,
#         ref="results/ref/{ref}/{ref}.fa",
#         regions="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
#         sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites",
#         idx="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites.idx",
#     output:
#         bcf=temp(
#             "results/datasets/{dataset}/bcfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.GLonly.bcf"
#         ),
#         maf=temp(
#             "results/datasets/{dataset}/bcfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.GLonly.mafs.gz"
#         ),
#         arg="results/datasets/{dataset}/bcfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.GLonly.arg",
#     log:
#         "logs/{dataset}/angsd/doBcf1/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.GLonly.log",
#     benchmark:
#         "benchmarks/{dataset}/angsd/doBcf1/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.GLonly.log"
#     container:
#         angsd.angsd_container
#     wildcard_constraints:
#         population="|".join(["all"]),
#         dp=".{0}",
#     params:
#         gl_model=config["params"]["angsd"]["gl_model"],
#         extra=config["params"]["angsd"]["extra"],
#         mapQ=config["mapQ"],
#         baseQ=config["baseQ"],
#         snp_pval=config["params"]["angsd"]["snp_pval"],
#         minind=lambda w: math.ceil(len(angsd.get_samples_from_pop(w.population)) * config["analyses"]["dataset_missing_data"]),
#         out=lambda w, output: os.path.splitext(output.arg)[0],
#     threads: lambda wildcards, attempt: attempt
#     resources:
#         runtime=lambda wildcards, attempt: attempt * 720,
#     shell:
#         """
#         angsd -doBcf 1 -bam {input.bam} -GL {params.gl_model} -ref {input.ref} \
#             -doMajorMinor 1 -doMaf 1 -SNP_pval {params.snp_pval} -minInd {params.minind} \
#             -nThreads {threads} {params.extra} \
#             -minMapQ {params.mapQ} -minQ {params.baseQ} -sites {input.sites} \
#             -dopost 1 --ignore-RG 0 -dogeno 1 -docounts 1 -rf {input.regions} \
#             -setMinDepthInd 10 -out {params.out} &> {log}
#         """


# rule concat_bcf:
#     """
#     Combine chromosome group GL BCFs into a full genome BCF.
#     """
#     input:
#         calls=lambda w: expand(
#             "results/datasets/{{dataset}}/bcfs/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.GLonly.bcf",
#             chunk=angsd.chunklist,
#         ),
#     output:
#         bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.GLonly.bcf",
#         pos="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.GLonly.pos",
#     log:
#         "logs/{dataset}/bcftools/concat/{dataset}.{ref}_{population}{dp}_{sites}-filts_merge-bcf.GLonly.log",
#     benchmark:
#         "benchmarks/{dataset}/bcftools/concat/{dataset}.{ref}_{population}{dp}_{sites}-filts_merge-bcf.GLonly.log"
#     conda:
#         "../envs/bcftools.yaml"
#     wildcard_constraints:
#         population="|".join(["all"]),
#         dp=".{0}",
#     resources:
#         runtime=lambda wildcards, attempt: attempt * 60,
#     shell:
#         """
#         bcftools concat -Ob -o {output.bcf} {input} 2> {log}
#         bcftools view {output.bcf} | grep -v "#" | cut -f1,2 > {output.pos} 2>> {log}
#         """


# rule bcftools_roh_likes:
#     """
#     Estimate runs of homozygosity with BCFtools using the genotype likelihoods
#     as input. Use allele frequencies from ANGSD (in BCF).
#     """
#     input:
#         "results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.GLonly.bcf",
#     output:
#         roh="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{population}{dp}_{sites}-filts.GLonly.roh",
#         sites="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{population}{dp}_{sites}-filts.sites.GLonly.roh",
#         regs="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{population}{dp}_{sites}-filts.regs.GLonly.roh",
#     log:
#         "logs/{dataset}/bcftools/roh/{dataset}.{ref}_{population}{dp}_{sites}-filts_roh.GLonly.log",
#     wildcard_constraints:
#         population="|".join(["all"]),
#         dp=".{0}",
#     conda:
#         "../envs/bcftools.yaml"
#     threads: lambda wildcards, attempt: attempt
#     params:
#         recrate=config["recrate"],
#     shell:
#         """
#         (bcftools roh -M {params.recrate} --AF-dflt 0.4 -o {output.roh} {input}
#         awk '$1=="ST"' {output.roh} > {output.sites}
#         awk '$1=="RG"' {output.roh} > {output.regs}) 2> {log}
#         """


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

#         while read chr; do
#             awk -v chr=$chr '$1 == chr' {input.bed} >> {output.bed}
#         done < {input.reg}
#         """


# rule bcftools_mpileup_call:
#     """
#     Call variants per population, calling genotypes.
#     """
#     input:
#         alignments=lambda w: expand(
#             "results/datasets/{{dataset}}/bams/{sample}.{{ref}}.bam",
#             sample=angsd.get_samples_from_pop(w.population),
#         ),
#         bai=lambda w: expand(
#             "results/datasets/{{dataset}}/bams/{sample}.{{ref}}.bam.bai",
#             sample=angsd.get_samples_from_pop(w.population),
#         ),
#         ref="results/ref/{ref}/{ref}.fa",
#         index="results/ref/{ref}/{ref}.fa.fai",
#         regions="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_chunk{chunk}_{sites}-filts.bed",
#     output:
#         bcf="results/datasets/{dataset}/bcfs/chunks/{dataset}.{ref}_{population}_chunk{chunk}_{sites}-filts.bcf",
#     conda:
#         "../envs/bcftools.yaml"
#     params:
#         extra=lambda w, input: f"--no-BAQ --min-MQ 30 --min-BQ 30 -a FORMAT/AD,FORMAT/DP,INFO/AD",
#     log:
#         "logs/{dataset}/bcftools/mpileup/{dataset}.{ref}_{population}_chunk{chunk}_{sites}-filts.log",
#     resources:
#         runtime=360,
#     shell:
#         """
#         bcftools mpileup --threads {threads} -f {input.ref} -R {input.regions} -Ou \
#             {params.extra} {input.alignments} |
#         bcftools call -m --variants-only --skip-variants indels -f GQ,GP -o {output.bcf}
#         """
# rule bcf_filter:
#     """
#     Filter called genotypes to exclude genotypes with < 30 genotype quality and
#     depth < 10.
#     """
#     input:
#         "results/datasets/{dataset}/bcfs/chunks/{dataset}.{ref}_{population}_chunk{chunk}_{sites}-filts.bcf",
#     output:
#         "results/datasets/{dataset}/vcfs/chunks/{dataset}.{ref}_{population}_chunk{chunk}_{sites}-filts.filter.vcf.gz",
#     conda:
#         "../envs/bcftools.yaml"
#     log:
#         "logs/{dataset}/bcftools/filter/{dataset}.{ref}_{population}_chunk{chunk}_{sites}-filts.log",
#     params:
#         miss= 1 - config["analyses"]["dataset_missing_data"]
#     shell:
#         """
#         bcftools filter -e'QUAL<30' -Ou {input} | bcftools +setGT -Ou -- -t q -n . \
#             -i'FMT/DP<5 | FMT/GQ<30' | bcftools +fill-tags -Ou -- -t all | \
#             bcftools filter -i 'F_MISSING<={params.miss}' -Oz -o {output} 2> {log}
#         """
# rule bcf_concat:
#     """
#     Combine per chromosome group BCFs for called genotypes into a full genome
#     VCF for each population.
#     """
#     input:
#         expand(
#             "results/datasets/{{dataset}}/vcfs/chunks/{{dataset}}.{{ref}}_{{population}}_chunk{chunk}_{{sites}}-filts.filter.vcf.gz",
#             chunk=angsd.chunklist,
#         ),
#     output:
#         vcf="results/datasets/{dataset}/vcfs/{dataset}.{ref}_{population}_{sites}-filts.filter.vcf.gz",
#         pos="results/datasets/{dataset}/vcfs/{dataset}.{ref}_{population}_{sites}-filts.filter.pos",
#     conda:
#         "../envs/bcftools.yaml"
#     log:
#         "logs/{dataset}/bcftools/concat/{dataset}.{ref}_{population}_{sites}-filts.log",
#     shell:
#         """
#         bcftools concat -Oz -o {output.vcf} {input} 2> {log}
#         zcat {output.vcf} | grep -v "#" | cut -f1,2 > {output.pos} 2>> {log}
#         """
# rule bcftools_roh_calls:
#     """
#     Estimate runs of homozygosity with BCFtools using the genotype calls as
#     input. Calculate allele-frequencies from genotype calls.
#     """
#     input:
#         "results/datasets/{dataset}/vcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filter.vcf.gz",
#     output:
#         roh="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{population}{dp}_{sites}-filts.calls.roh",
#         sites="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{population}{dp}_{sites}-filts.sites.calls.roh",
#         regs="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{population}{dp}_{sites}-filts.regs.calls.roh",
#     log:
#         "logs/{dataset}/bcftools/roh/{dataset}.{ref}_{population}{dp}_{sites}-filts_roh.calls.log",
#     wildcard_constraints:
#         population="|".join(["all"]),
#         dp=".{0}",
#     conda:
#         "../envs/bcftools.yaml"
#     threads: lambda wildcards, attempt: attempt
#     params:
#         recrate=config["recrate"],
#     shell:
#         """
#         (bcftools roh -G30 -M {params.recrate} --AF-dflt 0.4 -o {output.roh} {input}
#         awk '$1=="ST"' {output.roh} > {output.sites}
#         awk '$1=="RG"' {output.roh} > {output.regs}) 2> {log}
#         """
