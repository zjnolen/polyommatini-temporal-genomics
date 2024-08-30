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
        bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}.bamlist",
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
        minind=angsd.get_minind,
        mininddp=config["params"]["angsd"]["mindepthind"],
        out=lambda w, output: os.path.splitext(output.thetas)[0],
    threads: lambda w, attempt: attempt * 2
    resources:
        runtime=lambda w, attempt: attempt * 1200,
    shell:
        """
        (angsd -doSaf 1 -bam {input.sublist} -GL {params.gl_model} \
            -ref {input.ref} -nThreads {threads} {params.extra} \
            {params.minind} -minMapQ {params.mapQ} -minQ {params.baseQ} \
            -sites {input.sites} -anc {input.ref} -noTrans {params.trans} \
            -rf <(cut -f1 {input.sites} | sort | uniq | sed -e 's/$/:/') \
            -{params.extra_saf} -setMinDepthInd {params.mininddp} \
            -out {resources.tmpdir}/saf

        realSFS {resources.tmpdir}/saf.saf.idx -fold {params.fold} \
            -P {threads} > {output.sfs}

        realSFS saf2theta {resources.tmpdir}/saf.saf.idx -sfs {output.sfs} \
            -fold {params.fold} -outname {resources.tmpdir}/thetas
        
        thetaStat do_stat {resources.tmpdir}/thetas.thetas.idx \
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
            population=angsd.pop_list,
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
        (printf "pop\tdownsample.n\tdownsample.rep\tpi\twatterson\ttajima\n" > {output}
        cat {input} >> {output}) 2> {log}
        """
