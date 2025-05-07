rule bcftools_subset:
    input:
        vcf=config["gone_input_bcf"],
        tbi=config["gone_input_bcf"] + ".tbi",
    output:
        vcf="results/datasets/{dataset}/vcfs/{dataset}.{ref}_{population}_gone.vcf.gz",
        tbi="results/datasets/{dataset}/vcfs/{dataset}.{ref}_{population}_gone.vcf.gz.tbi",
    wildcard_constraints:
        population="|".join(angsd.pop_list),
    conda:
        "../envs/bcftools121.yaml"
    threads: 2
    resources:
        runtime="8h",
    params:
        samples=lambda w: ",".join(
            angsd.samples.index[
                angsd.samples.population == w.population
            ].values.tolist()
        ),
    shell:
        """
        bcftools view -Oz -s {params.samples} --force-samples {input.vcf} \
            > {output.vcf}
        tabix {output.vcf}
        """


rule bcf2ped:
    """
    Convert BCF to PED, the input for GONE
    """
    input:
        vcf="results/datasets/{dataset}/vcfs/{dataset}.{ref}_{population}_gone.vcf.gz",
        tbi="results/datasets/{dataset}/vcfs/{dataset}.{ref}_{population}_gone.vcf.gz.tbi",
    output:
        ped="results/datasets/{dataset}/peds/{dataset}.{ref}_{population}_gone.ped",
        map="results/datasets/{dataset}/peds/{dataset}.{ref}_{population}_gone.map",
        nosex=temp(
            "results/datasets/{dataset}/peds/{dataset}.{ref}_{population}_gone.nosex"
        ),
    container:
        "docker://quay.io/biocontainers/plink:1.90b6.18--h779adbc_1"
    params:
        pre=lambda w, output: os.path.splitext(output.ped)[0],
    resources:
        runtime=lambda wildcards, attempt: attempt * 360,
    shell:
        """
        plink --vcf {input.vcf} --allow-extra-chr --recode --out {params.pre}
        """


rule GONE:
    """
    Use population PED files to run GONE demographic analyses.
    """
    input:
        settings=config["GONE_input_params"],
        ped="results/datasets/{dataset}/peds/{dataset}.{ref}_{population}_gone.ped",
        map="results/datasets/{dataset}/peds/{dataset}.{ref}_{population}_gone.map",
    output:
        direct=directory(
            "results/datasets/{dataset}/analyses/gone/{dataset}.{ref}_{population}"
        ),
        ped="results/datasets/{dataset}/analyses/gone/{dataset}.{ref}_{population}/{dataset}.{ref}_{population}.ped",
        map="results/datasets/{dataset}/analyses/gone/{dataset}.{ref}_{population}/{dataset}.{ref}_{population}.map",
        d2="results/datasets/{dataset}/analyses/gone/{dataset}.{ref}_{population}/Output_d2_{dataset}.{ref}_{population}",
        out="results/datasets/{dataset}/analyses/gone/{dataset}.{ref}_{population}/OUTPUT_{dataset}.{ref}_{population}",
        Ne="results/datasets/{dataset}/analyses/gone/{dataset}.{ref}_{population}/Output_Ne_{dataset}.{ref}_{population}",
        timefile="results/datasets/{dataset}/analyses/gone/{dataset}.{ref}_{population}/timefile",
        seedfile="results/datasets/{dataset}/analyses/gone/{dataset}.{ref}_{population}/seedfile",
        hwd="results/datasets/{dataset}/analyses/gone/{dataset}.{ref}_{population}/outfileHWD",
    container:
        "docker://zjnolen/gone:20231231-134996d"
    threads: 20
    resources:
        runtime="12h",
    shell:
        """
        ln -sr {input.ped} {output.ped}
        ln -sr {input.map} {output.map}
        cp {input.settings} {output.direct}/INPUT_PARAMETERS_FILE
        cd {output.direct}
        sed -i "s/threads=-99/threads={threads}/g" INPUT_PARAMETERS_FILE
        script_GONE.sh {wildcards.dataset}.{wildcards.ref}_{wildcards.population}
        """


rule compile_GONE:
    input:
        expand(
            "results/datasets/{{dataset}}/analyses/gone/{{dataset}}.{{ref}}_{population}",
            population=list(
                set(angsd.pop_list)
                & set(
                    [
                        "ESkane2020",
                        "ESkane2021",
                        "ESkane2022",
                        "NSmaland2021",
                        "Oland2020",
                        "Oland2021",
                        "SESkane2022",
                        "SWSkane2021",
                        "WSkane2020",
                        "WSkane2021",
                    ]
                )
            ),
        ),
    output:
        "results/datasets/{dataset}/analyses/gone/STATUS_{dataset}.{ref}_populations",
    shell:
        """
        echo "All pops done!" > {output}
        """
