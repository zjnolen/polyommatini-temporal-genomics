localrules:
    calc_heterozygosity,
    pixy_poplist,


rule calc_heterozygosity_gencalls:
    """
    Calculate heterozygosity from the stats file of a BCF. Hets / total GTs
    """
    input:
        stats="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcf.stats",
    output:
        het="results/datasets/{dataset}/analyses/heterozygosity/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcftools.het.tsv",
    container:
        angsd.shell_container
    shell:
        """
        echo "sample\thetperbp\thetperalt" > {output.het}
        grep PSC {input.stats} | grep -v "#" | \
            awk '{{print $3"\t"$6/($6+$5+$4)"\t"$6/(2*$5+$6)}}' >> {output.het}
        """


rule pixy_poplist:
    """
    Generate a list of samples in the two column format expected by Pixy.
    Exclude samples that were excluded from calling
    """
    input:
        "results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
    output:
        "results/datasets/{dataset}/poplists/{dataset}_all.indiv.pixy.list",
    container:
        angsd.shell_container
    params:
        exclude=config["calling_drop"],
    shell:
        r"""
        echo "{params.exclude}" | tr ' ' '\n' | grep -v -f - {input} | \
            cut -f1-2 | tail -n+2 > {output}
        """


rule pixy:
    """
    Calculates Fst, nucleotide diversity, Dxy, in windows from called genotypes.
    """
    input:
        vcf="results/datasets/{dataset}/vcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.vcf.gz",
        tbi="results/datasets/{dataset}/vcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.vcf.gz.tbi",
        pops="results/datasets/{dataset}/poplists/{dataset}_all.indiv.pixy.list",
    output:
        fold=directory(
            "results/datasets/{dataset}/analyses/pixy/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}",
        ),
        pi="results/datasets/{dataset}/analyses/pixy/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}/pixy_pi.txt",
        fst="results/datasets/{dataset}/analyses/pixy/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}/pixy_fst.txt",
        dxy="results/datasets/{dataset}/analyses/pixy/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}/pixy_dxy.txt",
    params:
        winsize=50000,
    conda:
        "../envs/pixy.yaml"
    threads: 10
    resources:
        runtime="6h",
    shell:
        """
        pixy --stats pi fst dxy --vcf {input.vcf} --populations {input.pops} \
            --window_size {params.winsize} --fst_type hudson \
            --output_folder {output.fold}
        """
