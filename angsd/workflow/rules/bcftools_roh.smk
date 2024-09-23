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
    wildcard_constraints:
        miss=config["bcf_missing"],
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
    wildcard_constraints:
        miss=config["bcf_missing"],
    params:
        bins=config["params"]["ngsf-hmm"]["roh_bins"],
        minroh=config["params"]["ngsf-hmm"]["min_roh_length"],
        outpre=lambda w, output: output["barplot"].removesuffix(".froh_bins.svg"),
        roh_phred=config["min_phred_roh"],
    script:
        "../scripts/plot_Froh.R"


rule bcftools_roh_testAF:
    input:
        "results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.bcf",
    output:
        roh="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.AFtest{AF}.{homref}.roh",
        sites="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.AFtest{AF}.{homref}.sites.roh",
        regs="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.AFtest{AF}.{homref}.regs.roh",
    log:
        "logs/{dataset}/bcftools/roh/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.AFtest{AF}.{homref}.log",
    conda:
        "../envs/bcftools.yaml"
    wildcard_constraints:
        miss=config["bcf_missing"],
        homref="--ignore-homref|",
    threads: lambda wildcards, attempt: attempt
    params:
        recrate=config["recrate"],
    shell:
        """
        (bcftools roh -G30 {wildcards.homref} -M {params.recrate} \
            --AF-dflt {wildcards.AF} -o {output.roh} {input}
        awk '$1=="ST"' {output.roh} > {output.sites}
        awk '$1=="RG"' {output.roh} > {output.regs}) 2> {log}
        """


rule bcftools_roh_testAF_plot:
    input:
        roh="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.AFtest{AF}.{homref}.regs.roh",
        inds="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
        autos=angsd.get_auto_sum,
    output:
        barplot="results/datasets/{dataset}/plots/inbreeding/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.AFtest{AF}.{homref}.bcftools.froh_bins.svg",
        scatter="results/datasets/{dataset}/plots/inbreeding/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.AFtest{AF}.{homref}.bcftools.cumroh_nroh.svg",
        froh="results/datasets/{dataset}/plots/inbreeding/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.AFtest{AF}.{homref}.bcftools.ind_froh.tsv",
    log:
        "logs/{dataset}/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts_filtered_mindp{mindp}-biallelic.{trans}.fmiss{miss}.AFtest{AF}.{homref}_plot.log",
    conda:
        "../envs/r.yaml"
    wildcard_constraints:
        miss=config["bcf_missing"],
        homref="--ignore-homref|",
    params:
        bins=config["params"]["ngsf-hmm"]["roh_bins"],
        minroh=config["params"]["ngsf-hmm"]["min_roh_length"],
        outpre=lambda w, output: output["barplot"].removesuffix(".froh_bins.svg"),
        roh_phred=config["min_phred_roh"],
    script:
        "../scripts/plot_Froh.R"
