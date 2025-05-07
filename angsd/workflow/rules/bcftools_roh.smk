rule bcftools_roh:
    """
    Estimate runs of homozygosity from called genotypes using bcftools.
    """
    input:
        "results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcf",
    output:
        roh="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.roh",
        sites="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.sites.roh",
        regs="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.regs.roh",
    log:
        "logs/{dataset}/bcftools/roh/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.log",
    conda:
        "../envs/bcftools121.yaml"
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
    """
    Plot the results of the bcftools ROH analysis for visualization.
    """
    input:
        roh="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.regs.roh",
        inds="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
        autos=angsd.get_auto_sum,
    output:
        barplot=report(
            "results/datasets/{dataset}/plots/inbreeding/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcftools.froh_bins.svg",
            category="06 Inbreeding",
            labels=lambda w: {
                "Filter": "{sites}",
                **angsd.dp_report(w),
                "minCallDP": "{mindp}",
                "Type": "Froh Bins Barplot (bcftools)",
            },
        ),
        scatter=report(
            "results/datasets/{dataset}/plots/inbreeding/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcftools.cumroh_nroh.svg",
            category="06 Inbreeding",
            labels=lambda w: {
                "Filter": "{sites}",
                **angsd.dp_report(w),
                "minCallDP": "{mindp}",
                "Type": "Nroh ~ Lroh Scatterplot (bcftools)",
            },
        ),
        froh="results/datasets/{dataset}/plots/inbreeding/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcftools.ind_froh.tsv",
    log:
        "logs/{dataset}/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts_filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}_plot.log",
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
