# These rules downsample the various diversity statistics used to lower sample
# sizes to assess how small a sample size can approach the full sample size.

# adapted from original version at workflow/rules/2.0_downsampling.smk in
# https://zenodo.org/records/13383389


localrules:
    downsample_bamlist,


rule downsample_bamlist:
    """
    Generate downsampled lists of bam files, these will be used as input to
    ANGSD to calculate diversity statistics for downsampled replicates.
    """
    input:
        bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=angsd.get_bamlist_bams,
        bais=angsd.get_bamlist_bais,
    output:
        sublist="results/datasets/{dataset}/bamlists/downsampled/{dataset}.{ref}_{population}{dp}.N{samplesize}-rep{rep}.bamlist",
    log:
        "logs/{dataset}/downsampled/sublist/{dataset}.{ref}_{population}{dp}.N{samplesize}-rep{rep}.log",
    benchmark:
        "benchmarks/{dataset}/downsampled/sublist/{dataset}.{ref}_{population}{dp}.N{samplesize}-rep{rep}.log"
    conda:
        "../envs/shell.yaml"
    resources:
        runtime="10m",
    shell:
        """
        shuf -n{wildcards.samplesize} {input.bamlist} > {output.sublist} 2> {log}
        """


rule downsample_thetas:
    """
    Estimate pi, theta, and Tajima's D using ANGSD on the downsampled bam lists.
    """
    input:
        sublist="results/datasets/{dataset}/bamlists/downsampled/{dataset}.{ref}_{population}{dp}.N{samplesize}-rep{rep}_{sites}-filts.bamlist",
        ref="results/ref/{ref}/{ref}.fa",
        sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites",
        idx="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites.idx",
    output:
        sfs=temp(
            ensure(
                "results/datasets/{dataset}/analyses/sfs/downsampled/{dataset}.{ref}_{population}{dp}.N{samplesize}-rep{rep}_{sites}-filts.thetaWindows.{win}_{step}.sfs",
                non_empty=True,
            )
        ),
        thetas=temp(
            "results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population}{dp}.N{samplesize}-rep{rep}_{sites}-filts.thetaWindows.{win}_{step}.pestPG"
        ),
    container:
        angsd.angsd_container
    log:
        "logs/{dataset}/downsampled/thetas/{dataset}.{ref}_{population}{dp}.N{samplesize}-rep{rep}_{sites}-filts.thetaWindows.{win}_{step}.log",
    group:
        "theta"
    benchmark:
        "benchmarks/{dataset}/downsampled/thetas/{dataset}.{ref}_{population}{dp}.N{samplesize}-rep{rep}_{sites}-filts.thetaWindows.{win}_{step}.log"
    params:
        gl_model=config["params"]["angsd"]["gl_model"],
        extra=config["params"]["angsd"]["extra"],
        extra_saf=config["params"]["angsd"]["extra_saf"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        trans=angsd.get_trans,
        minind=lambda w: int(
            int(w.samplesize) * config["params"]["angsd"]["minind_dataset"]
        ),
        mininddp=config["params"]["angsd"]["mindepthind"],
        fold=config["params"]["realsfs"]["fold"],
        pre=lambda w: f"{w.dataset}.{w.ref}_{w.population}{w.dp}.N{w.samplesize}-rep{w.rep}_{w.sites}-filts",
        out=lambda w, output: os.path.splitext(output.thetas)[0],
    threads: lambda w, attempt: attempt * 2
    resources:
        runtime=lambda w, attempt: attempt * 1200,
    shell:
        """
        (angsd -doSaf 1 -bam {input.sublist} -GL {params.gl_model} \
            -ref {input.ref} -nThreads {threads} {params.extra} \
            -minInd {params.minind} -minMapQ {params.mapQ} \
            -minQ {params.baseQ} -sites {input.sites} -anc {input.ref} \
            -noTrans {params.trans} \
            -rf <(cut -f1 {input.sites} | uniq | sed -e 's/$/:/') \
            {params.extra_saf} -setMinDepthInd {params.mininddp} \
            -out {resources.tmpdir}/{params.pre}

        realSFS {resources.tmpdir}/{params.pre}.saf.idx -fold {params.fold} \
            -P {threads} > {output.sfs}

        realSFS saf2theta {resources.tmpdir}/{params.pre}.saf.idx -sfs {output.sfs} \
            -fold {params.fold} -outname {resources.tmpdir}/{params.pre}
        
        thetaStat do_stat {resources.tmpdir}/{params.pre}.thetas.idx \
            -win {wildcards.win} -type 0 -step {wildcards.step} \
            -outnames {params.out}) &> {log}
        """


rule average_downsampled_thetas:
    """
    Get a mean estimate for each downsampled pi, theta, and Tajima's D from the
    sliding window estimates.
    """
    input:
        "results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population}{dp}.N{samplesize}-rep{rep}_{sites}-filts.thetaWindows.{win}_{step}.pestPG",
    output:
        ensure(
            "results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population}{dp}.N{samplesize}-rep{rep}_{sites}-filts.thetaMean.{win}_{step}.tsv"
        ),
        "results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population}{dp}.N{samplesize}-rep{rep}_{sites}-filts.windows.{win}_{step}.tsv",
    conda:
        "../envs/r.yaml"
    group:
        "theta"
    log:
        "logs/{dataset}/thetaStat/average/downsampled/{dataset}.{ref}_{population}{dp}.N{samplesize}-rep{rep}_{sites}-filts.thetaMean.{win}_{step}.log",
    params:
        minsites=config["params"]["thetas"]["minsites"],
    script:
        "../scripts/average_downsampled_thetas.R"


rule aggregate_downsampled_thetas:
    """
    Compile downsampled mean pi, theta and Tajima's D into a table with all
    downsampled replicates.
    """
    input:
        expand(
            "results/datasets/{{dataset}}/analyses/thetas/downsampled/{{dataset}}.{{ref}}_{population}{{dp}}.N{samplesize}-rep{rep}_{{sites}}-filts.thetaMean.{{win}}_{{step}}.tsv",
            population=[
                pop
                for pop in angsd.pop_list
                if pop not in ["ESkane1917", "WSkane1943"]
            ],
            samplesize=config["downsample_sizes"],
            rep=list(range(0, config["downsample_reps"][0])),
        ),
    output:
        "results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_all{dp}.downsampled_{sites}-filts.thetaMean.{win}_{step}.tsv",
    log:
        "logs/{dataset}/thetaStat/aggregate/downsampled/{dataset}.{ref}_all{dp}.downsampled_{sites}-filts.thetaMean.{win}_{step}.log",
    conda:
        "../envs/shell.yaml"
    resources:
        runtime=lambda w, attempt: attempt * 360,
    shell:
        """
        (printf "pop\tdownsample.n\tdownsample.rep\tpi\tpi.lower\tpi.upper\twatterson\twatterson.lower\twatterson.upper\ttajima\ttajima.lower\ttajima.upper\n" > {output}
        cat {input} >> {output}) 2> {log}
        """


rule wilcoxon_test_thetas:
    """
    Performs a paired Wilcoxon ranked sign test on theta distributions between
    populations after subsampling individuals. For each subsampling replicate it
    reports the median difference between modern and historical thetas across
    windows, and the p-value, Z-score, and effect size of the test.
    """
    input:
        pop1=expand(
            "results/datasets/{{dataset}}/analyses/thetas/downsampled/{{dataset}}.{{ref}}_{{population1}}{{dp}}.N{{samplesize}}-rep{rep}_{{sites}}-filts.windows.{win}_{step}.tsv",
            rep=list(range(0, config["downsample_reps"][0])),
            win=config["params"]["thetas"]["win_size"],
            step=config["params"]["thetas"]["win_step"],
        ),
        pop2=expand(
            "results/datasets/{{dataset}}/analyses/thetas/downsampled/{{dataset}}.{{ref}}_{{population2}}{{dp}}.N{{samplesize}}-rep{rep}_{{sites}}-filts.windows.{win}_{step}.tsv",
            rep=list(range(0, config["downsample_reps"][0])),
            win=config["params"]["thetas"]["win_size"],
            step=config["params"]["thetas"]["win_step"],
        ),
    output:
        stats="results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population1}-{population2}{dp}.N{samplesize}_{sites}-filts.wilcoxon.tsv",
    log:
        "logs/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population1}-{population2}{dp}.N{samplesize}_{sites}-filts.wilcoxon.log",
    conda:
        "../envs/r-wilcoxon.yaml"
    params:
        nreps=config["downsample_reps"],
    script:
        "../scripts/theta_wilcoxon.R"


def get_temp_wilcox(wildcards):
    samples = angsd.samples
    samples = samples[samples["population"] != "ESkane1917"]
    samples = samples[samples["population"] != "WSkane1943"]
    samples["year"] = samples["population"].str.slice(start=-4)
    samples["region"] = samples["population"].str.slice(stop=-4)
    hisset = set(samples[samples["time"] == "historical"]["region"])
    modset = set(samples[samples["time"] == "modern"]["region"])
    regions = hisset.intersection(modset)
    modpop = sorted(
        list(
            set(
                samples[
                    (samples["time"] == "modern") & (samples["region"].isin(regions))
                ]["population"]
            )
        )
    )
    hispop = sorted(
        list(
            set(
                samples[
                    (samples["time"] == "historical")
                    & (samples["region"].isin(regions))
                ]["population"]
            )
        )
    )
    if "SESkane" in modset:
        modpop = modpop + ["SESkane2022"]
        hispop = hispop + ["ESkane1956"]
    return expand(
        "results/datasets/{{dataset}}/analyses/thetas/downsampled/{{dataset}}.{{ref}}_{population1}-{population2}{{dp}}.N{{samplesize}}_{{sites}}-filts.wilcoxon.tsv",
        zip,
        population1=modpop,
        population2=hispop,
    )


rule merge_temporal_wilcoxon:
    """
    Creates a data table of all Wilcoxon ranked sign test outputs, grabbing only
    direct comparions of historical and modern populations so that differences
    in the table are modern - historical.
    """
    input:
        get_temp_wilcox,
    output:
        "results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_modern-historical{dp}.N{samplesize}_{sites}-filts.wilcoxon.tsv",
    conda:
        "../envs/shell.yaml"
    shell:
        """
        printf "modpop\thistpop\tdownsample.n\tdownsample.rep\tpi.diff\tpi.retained\tpi.z\tpi.effsize\tpi.pval\twatterson.diff\twatterson.retained\twatterson.z\twatterson.effsize\twatterson.pval\ttajima.diff\ttajima.retained\ttajima.z\ttajima.effsize\ttajima.pval\n" > {output}
        cat {input} >> {output}
        """
