rule angsd_doIBS_refmajor:
    """
    Generates an identity by state file with the reference allele as the major
    allele. This randomly samples a base for each individual and assigns 0 if
    it is major and 1 if it is not. This can then be used to calculate an IBS
    distance from the reference.
    """
    input:
        bam="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=angsd.get_bamlist_bams,
        bais=angsd.get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
        reffai="results/ref/{ref}/{ref}.fa.fai",
    output:
        ibs=temp("results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_rmtrans{rmtrans}.ibs.gz"),
        arg="results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_rmtrans{rmtrans}.arg",
        stats="results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_rmtrans{rmtrans}.refibs.tsv"
    wildcard_constraints:
        population="|".join(angsd.samples.index),
    log:
        "logs/{dataset}/angsd/doIBS_refmajor/{dataset}.{ref}_{population}{dp}_rmtrans{rmtrans}.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doIBS_refmajor/{dataset}.{ref}_{population}{dp}_rmtrans{rmtrans}.log"
    container:
        angsd.angsd_container
    params:
        extra=config["params"]["angsd"]["extra"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        r"""
        angsd -doIBS 1 -bam {input.bam} -ref {input.ref} -nThreads {threads} \
            -doMajorMinor 4 {params.extra} -GL 2 -minMapQ {params.mapQ} \
            -minQ {params.baseQ} -doCounts 1 -output01 1 \
            -out {params.out} -rmtrans {wildcards.rmtrans} &> {log}
        
        ibs=$(zcat {output.ibs} | tail -n+2 | \
            awk '{{ sum += $5 }} END {{ if (NR > 0) print 1 - (sum / NR) }}')
        echo "{wildcards.population}\t$ibs" > {output.stats} 2> {log}
        """


rule angsd_doIBS_refmajor_sites:
    """
    Generates an identity by state file with the reference allele as the major
    allele. This randomly samples a base for each individual and assigns 0 if
    it is major and 1 if it is not. This can then be used to calculate an IBS
    distance from the reference. Only does it on filtered sites.
    """
    input:
        sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites",
        idx="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites.idx",
        bam="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=angsd.get_bamlist_bams,
        bais=angsd.get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
        reffai="results/ref/{ref}/{ref}.fa.fai",
    output:
        ibs=temp("results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_{sites}-filts_rmtrans{rmtrans}.ibs.gz"),
        arg="results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_{sites}-filts_rmtrans{rmtrans}.arg",
        stats="results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_{sites}-filts_rmtrans{rmtrans}.refibs.tsv"
    wildcard_constraints:
        population="|".join(angsd.samples.index),
    log:
        "logs/{dataset}/angsd/doIBS_refmajor/{dataset}.{ref}_{population}{dp}_{sites}-filts_rmtrans{rmtrans}.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doIBS_refmajor/{dataset}.{ref}_{population}{dp}_{sites}-filts_rmtrans{rmtrans}.log"
    container:
        angsd.angsd_container
    params:
        extra=config["params"]["angsd"]["extra"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        r"""
        angsd -doIBS 1 -bam {input.bam} -ref {input.ref} -nThreads {threads} \
            -doMajorMinor 4 {params.extra} -GL 2 -minMapQ {params.mapQ} \
            -minQ {params.baseQ} -doCounts 1 -output01 1 -sites {input.sites} \
            -out {params.out} -rmtrans {wildcards.rmtrans} &> {log}
        
        ibs=$(zcat {output.ibs} | tail -n+2 | \
            awk '{{ sum += $5 }} END {{ if (NR > 0) print 1 - (sum / NR) }}')
        echo "{wildcards.population}\t$ibs" > {output.stats} 2> {log}
        """


rule merge_doIBS_refmajor:
    """
    Combines IBS files across chunks
    """
    input:
        lambda w: expand(
            "results/datasets/{{dataset}}/qc/ibs_refbias/{{dataset}}.{{ref}}_{population}{{dp}}_rmtrans{{rmtrans}}.refibs.tsv",
            population=angsd.samples.index,
        ),
    output:
        ibs="results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_all{dp}_rmtrans{rmtrans}.refibs.tsv",
    conda:
        "../envs/shell.yaml"
    resources:
        runtime=lambda wildcards, attempt: attempt * 60,
    shell:
        r"""
        cat {input} > {output}
        """


rule merge_doIBS_refmajor_sites:
    """
    Combines IBS files across chunks
    """
    input:
        lambda w: expand(
            "results/datasets/{{dataset}}/qc/ibs_refbias/{{dataset}}.{{ref}}_{population}{{dp}}_{{sites}}-filts_rmtrans{{rmtrans}}.refibs.tsv",
            population=angsd.samples.index,
        ),
    output:
        ibs="results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_all{dp}_{sites}-filts_rmtrans{rmtrans}.refibs.tsv",
    conda:
        "../envs/shell.yaml"
    resources:
        runtime=lambda wildcards, attempt: attempt * 60,
    shell:
        r"""
        cat {input} > {output}
        """