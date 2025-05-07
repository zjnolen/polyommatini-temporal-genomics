"""
This Snakefile contains all the rules related to constructing the variant graph
with vg and then mapping the samples to it. These rules are run only
when the variant graph is constructed, which is a single Snakemake run using
only the modern samples for a species, which have been mapped with bwa mem. This
is done with a separate snakefile `Snakefile_vg-construct` to the main workflow.
"""


rule bcftools_joint_call4vg_chunk:
    """
    Calls genotypes for use in constructing a species-specific variation graph.

    Calls genotypes jointly across all samples in a species dataset, excluding
    those that have been dropped from genotype call analyses. Uses bcftools
    multiallelic caller and groups individuals by sample population for the
    calling model's HWE assumption. After calling, drop low quality positions
    (QUAL < 30) and genotypes (GQ < 30, DP < 6). Filters to only biallelic SNPs
    with a MAF >= 0.05 and < 0.4 missing data. This ensures that only relatively
    common SNPs are used in the variation graph.
    """
    input:
        bams=expand(
            "results/datasets/{{dataset}}/bams/{population}.{{ref}}{{dp}}.bam",
            population=[
                sample
                for sample in angsd.samples.index.tolist()
                if sample not in config["calling_drop"]
            ],
        ),
        bais=expand(
            "results/datasets/{{dataset}}/bams/{population}.{{ref}}{{dp}}.bam.bai",
            population=[
                sample
                for sample in angsd.samples.index.tolist()
                if sample not in config["calling_drop"]
            ],
        ),
        ref="results/ref/{ref}/{ref}.fa",
        index="results/ref/{ref}/{ref}.fa.fai",
        regions="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
        sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.bed",
        poplist="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
    output:
        bcf=temp(
            "results/datasets/{dataset}/bcfs/chunks/{chunk}/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-jointcall.bcf"
        ),
        idx=temp(
            "results/datasets/{dataset}/bcfs/chunks/{chunk}/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-jointcall.bcf.csi"
        ),
    conda:
        "../envs/bcftools121.yaml"
    threads: 2
    resources:
        runtime="8h",
    shell:
        """
        bcftools mpileup -f {input.ref} -R {input.regions} -T {input.sites} \
            -a "FORMAT/QS,FORMAT/AD,FORMAT/DP,INFO/AD" -B --min-MQ 30 \
            --min-BQ 30 -Ou {input.bams} | \
            bcftools call -m -a GQ,GP -G {input.poplist} -Ou | \
            bcftools filter -g 5 -i'QUAL >= 30' -Ou | \
            bcftools filter -i'FMT/GQ >= 30' -S . -Ou | \
            bcftools filter -i'FMT/DP >= {wildcards.mindp}' -S . -Ou | \
            bcftools view -V indels -M2 -Ou | \
            bcftools +fill-tags -Ou -- -t all | \
            bcftools filter -e'MAF<0.05 | F_MISSING>0.4' -Ob > {output.bcf}
        bcftools index -o {output.idx} {output.bcf}
        """


rule bcftools_concat_vg_chunks:
    """
    Concatenates BCF chunks into a single BCF with whole genome.
    """
    input:
        bcfs=expand(
            "results/datasets/{{dataset}}/bcfs/chunks/{chunk}/{{dataset}}.{{ref}}_all_allsites-filts.filtered_mindp6-allsites-jointcall.bcf",
            chunk=angsd.chunklist,
        ),
        idx=expand(
            "results/datasets/{{dataset}}/bcfs/chunks/{chunk}/{{dataset}}.{{ref}}_all_allsites-filts.filtered_mindp6-allsites-jointcall.bcf.csi",
            chunk=angsd.chunklist,
        ),
    output:
        vcf="results/datasets/{dataset}/vcfs/{dataset}.{ref}_all_allsites-filts.vcf.gz",
        idx="results/datasets/{dataset}/vcfs/{dataset}.{ref}_all_allsites-filts.vcf.gz.tbi",
        stats="results/datasets/{dataset}/vcfs/{dataset}.{ref}_all_allsites-filts.vcf.stats",
    conda:
        "../envs/bcftools121.yaml"
    threads: 4
    resources:
        runtime="6h",
    shell:
        """
        bcftools concat --threads {threads} -Oz {input.bcfs} > {output.vcf}
        tabix {output.vcf}
        bcftools stats -s - {output.vcf} > {output.stats}
        """


if len(angsd.samples.index.tolist()) > 0:

    rule construct_vg:
        """
        Construct a variant graph from the reference genome and the variants
        found in the modern samples.
        """
        input:
            ref="results/ref/{ref}/{ref}.fa",
            fai="results/ref/{ref}/{ref}.fa.fai",
            vcf=expand(
                "results/datasets/{dataset}/vcfs/{dataset}.{{ref}}_all_allsites-filts.vcf.gz",
                dataset=config["dataset"],
            ),
            idx=expand(
                "results/datasets/{dataset}/vcfs/{dataset}.{{ref}}_all_allsites-filts.vcf.gz.tbi",
                dataset=config["dataset"],
            ),
        output:
            vg="results/ref/{ref}/{ref}.modvars.vg",
        conda:
            "../envs/vg.yaml"
        threads: 8
        resources:
            runtime="12h",
        shell:
            """
            vg construct -r {input.ref} -v {input.vcf} -p > {output.vg}
            """

else:

    rule construct_vg:
        """
        Construct a variant graph from the reference genome. This version is
        for when there are no modern samples to call genotypes from (i.e., for
        Pl. argyrognomon and Sc. orion).
        """
        input:
            ref="results/ref/{ref}/{ref}.fa",
            fai="results/ref/{ref}/{ref}.fa.fai",
        output:
            vg="results/ref/{ref}/{ref}.modvars.vg",
        conda:
            "../envs/vg.yaml"
        threads: 8
        resources:
            runtime="12h",
        shell:
            """
            vg construct -r {input.ref} -p > {output.vg}
            """


rule index_vg:
    """
    Indexes the variant graph.
    """
    input:
        vg="results/ref/{ref}/{ref}.modvars.vg",
    output:
        prune="results/ref/{ref}/{ref}.modvars.prune.vg",
        xg="results/ref/{ref}/{ref}.modvars.xg",
        gcsa="results/ref/{ref}/{ref}.modvars.gcsa",
    conda:
        "../envs/vg.yaml"
    threads: 128
    resources:
        runtime="1d",
        slurm_partition="main",
        # mem_mb=450560,
    shell:
        """
        vg index -t {threads} -p -x {output.xg} {input.vg}
        vg prune -t {threads} -p -r {input.vg} > {output.prune}
        vg index -t {threads} -p -g {output.gcsa} {output.prune}
        """
