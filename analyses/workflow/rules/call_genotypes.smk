"""
Rules related to genotype calling in this project. These focus on using joint
calling to try to increase confidence in low depth genotypes, using missingness
filters independently for sample categories (historical/modern) to ensure that
missingness is even across treatments, and allelic balance based genotype
filters to reduce calls from rare mismatches that may be mapping errors.
Individual level calling is also possible, with rules towards the bottom.
"""


rule bcftools_joint_call_chunk:
    """
    Calls genotypes jointly across all samples in a species dataset for a region
    of the genome, excluding samples that have been dropped from genotype call
    analyses. Uses bcftools multiallelic caller and groups individuals by sample
    population for the calling model's HWE assumption. After calling, drop low
    quality positions (QUAL) and set genotypes with low depth and/or poor
    allelic balance to missing. Keeps monomorphic positions and biallelic SNPs
    and outputs BCFs both with and without transitions where DNA damage might be
    present. Applies no missing data filter yet, aside from how missingness
    informs QUAL score.
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
        bcf_trans=temp(
            "results/datasets/{dataset}/bcfs/chunks/{chunk}/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-jointcall_allbal{ablow}-{abhi}.trans.bcf"
        ),
        idx_trans=temp(
            "results/datasets/{dataset}/bcfs/chunks/{chunk}/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-jointcall_allbal{ablow}-{abhi}.trans.bcf.csi"
        ),
        bcf_notrans=temp(
            "results/datasets/{dataset}/bcfs/chunks/{chunk}/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-jointcall_allbal{ablow}-{abhi}.notrans.bcf"
        ),
        idx_notrans=temp(
            "results/datasets/{dataset}/bcfs/chunks/{chunk}/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-jointcall_allbal{ablow}-{abhi}.notrans.bcf.csi"
        ),
    conda:
        "../envs/bcftools121.yaml"
    threads: 2
    resources:
        runtime="8h",
    shell:
        """
        bcftools mpileup -f {input.ref} -R {input.regions} -T {input.sites} \
            -a "FORMAT/QS,FORMAT/AD,FORMAT/DP,INFO/AD" -B \
            --min-MQ 30 --min-BQ 20 -Ou {input.bams} | \
            bcftools call -m -a GQ,GP -G {input.poplist} -Ou | \
            bcftools filter -g 5 -i'QUAL >= 30' -Ou | \
            bcftools view -V indels -M2 -Ou | \
            bcftools +fill-tags -Ou -- -t all | \
            bcftools +setGT -Ou -- -t q -n . -i"FMT/DP<{wildcards.mindp}" | \
            bcftools +setGT -Ou -- -t q -n . -i'GT="het" & (FMT/VAF < {wildcards.ablow} | FMT/VAF > {wildcards.abhi})' | \
            bcftools +fill-tags -Ob -- -t all > {output.bcf_trans}
        bcftools view {output.bcf_trans} | \
            awk -F '\t' '!(($4 == "A" && $5 == "G") || ($4 == "G" && $5 == "A") || ($4 == "C" && $5 == "T") || ($4 == "T" && $5 == "C"))' | \
            bcftools view -Ob > {output.bcf_notrans}
        bcftools index -o {output.idx_trans} {output.bcf_trans}
        bcftools index -o {output.idx_notrans} {output.bcf_notrans}
        """


rule bcftools_concat_chunks:
    """
    Concatenates BCF chunks into a single BCF with whole genome.
    """
    input:
        bcfs=expand(
            "results/datasets/{{dataset}}/bcfs/chunks/{chunk}/{{dataset}}.{{ref}}_all{{dp}}_{{sites}}-filts.filtered_mindp{{mindp}}-allsites-{{call}}_allbal{{ablow}}-{{abhi}}.{{trans}}.bcf",
            chunk=angsd.chunklist,
        ),
        idx=expand(
            "results/datasets/{{dataset}}/bcfs/chunks/{chunk}/{{dataset}}.{{ref}}_all{{dp}}_{{sites}}-filts.filtered_mindp{{mindp}}-allsites-{{call}}_allbal{{ablow}}-{{abhi}}.{{trans}}.bcf.csi",
            chunk=angsd.chunklist,
        ),
    output:
        bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-{call}_allbal{ablow}-{abhi}.{trans}.bcf",
        idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-{call}_allbal{ablow}-{abhi}.{trans}.bcf.csi",
        stats="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-{call}_allbal{ablow}-{abhi}.{trans}.bcf.stats",
        psc="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-{call}_allbal{ablow}-{abhi}.{trans}.bcf.psc",
    conda:
        "../envs/bcftools121.yaml"
    threads: 4
    resources:
        runtime="6h",
    shell:
        """
        bcftools concat --threads {threads} -Ob {input.bcfs} > {output.bcf}
        bcftools index -o {output.idx} {output.bcf}
        bcftools stats -s - {output.bcf} > {output.stats}
        grep PSC {output.stats} > {output.psc}
        """


rule bcftools_biallelic_snps:
    """
    Filters a variant + invariant site bcf to biallelic snps, no MAF filter.
    """
    input:
        bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-{call}_allbal{ablow}-{abhi}.{trans}.bcf",
        idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-{call}_allbal{ablow}-{abhi}.{trans}.bcf.csi",
    output:
        bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.bcf",
        idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.bcf.csi",
        stats="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-biallelic-{call}_allbal{ablow}-{abhi}.{trans}.bcf.stats",
    wildcard_constraints:
        trans="notrans|trans",
    conda:
        "../envs/bcftools121.yaml"
    threads: 2
    resources:
        runtime="4h",
    shell:
        r"""
        bcftools view -v snps -m 2 -M 2 -i'MAF>0' -Ob {input.bcf} > {output.bcf}
        bcftools index -o {output.idx} {output.bcf}
        bcftools stats -s - {output.bcf} > {output.stats}
        """


if (
    len(angsd.samples.index[angsd.samples.time == "historical"].values.tolist()) > 0
    and len(angsd.samples.index[angsd.samples.time == "modern"].values.tolist()) > 0
):

    rule bcftools_missingness:
        """
        Filters based on a missingness threshold that must be met in both the
        historical and modern sample subsets.
        """
        input:
            bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-{gts}-{call}_allbal{ablow}-{abhi}.{trans}.bcf",
            idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-{gts}-{call}_allbal{ablow}-{abhi}.{trans}.bcf.csi",
        output:
            bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-{gts}-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcf",
            idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-{gts}-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcf.csi",
            stats="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-{gts}-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcf.stats",
            modsites=temp(
                "results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-{gts}-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcf.modsites"
            ),
            histsites=temp(
                "results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-{gts}-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcf.histsites"
            ),
        wildcard_constraints:
            trans="notrans|trans",
            gts="allsites|biallelic",
        conda:
            "../envs/bcftools121.yaml"
        threads: 6
        resources:
            runtime="6h",
        params:
            histsamps=",".join(
                angsd.samples.index[angsd.samples.time == "historical"].values.tolist()
            ),
            modsamps=",".join(
                angsd.samples.index[angsd.samples.time == "modern"].values.tolist()
            ),
        shell:
            """
            bcftools view -s {params.histsamps} --force-samples -Ou {input.bcf} | \
                bcftools +fill-tags -Ou -- -t F_MISSING | \
                bcftools query -i'F_MISSING < {wildcards.miss}' \
                    -f '%CHROM\t%POS' > {output.histsites}
            bcftools view -s {params.modsamps} --force-samples -Ou {input.bcf} | \
                bcftools +fill-tags -Ou -- -t F_MISSING | \
                bcftools query -i'F_MISSING < {wildcards.miss}' \
                    -f '%CHROM\t%POS' > {output.modsites}
            bcftools view -T {output.histsites} -Ou {input.bcf} | \
                bcftools view -T {output.modsites} -Ob > {output.bcf}
            bcftools stats -s - {output.bcf} > {output.stats}
            bcftools index -o {output.idx} {output.bcf}
            """

else:

    rule bcftools_missingness:
        """
        Filters based on a missingness threshold when BCF includes only modern
        or only historical samples.
        """
        input:
            bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-{gts}-{call}_allbal{ablow}-{abhi}.{trans}.bcf",
            idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-{gts}-{call}_allbal{ablow}-{abhi}.{trans}.bcf.csi",
        output:
            bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-{gts}-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcf",
            idx="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-{gts}-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcf.csi",
            stats="results/datasets/{dataset}/bcfs/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-{gts}-{call}_allbal{ablow}-{abhi}.{trans}.fmiss{miss}.bcf.stats",
        wildcard_constraints:
            trans="notrans|trans",
            gts="allsites|biallelic",
        conda:
            "../envs/bcftools121.yaml"
        threads: 6
        resources:
            runtime="6h",
        shell:
            """
            bcftools view -Ou {input.bcf} | \
                bcftools +fill-tags -Ou -- -t F_MISSING | \
                bcftools view -i'F_MISSING < {wildcards.miss}' -Ob \
                    > {output.bcf}
            bcftools stats -s - {output.bcf} > {output.stats}
            bcftools index -o {output.idx} {output.bcf}
            """


rule bcf2vcf:
    """
    Makes a bcf a vcf when needed.
    """
    input:
        bcf="results/datasets/{dataset}/bcfs/{prefix}.bcf",
        idx="results/datasets/{dataset}/bcfs/{prefix}.bcf.csi",
    output:
        vcf="results/datasets/{dataset}/vcfs/{prefix}.vcf.gz",
        tbi="results/datasets/{dataset}/vcfs/{prefix}.vcf.gz.tbi",
    conda:
        "../envs/bcftools121.yaml"
    threads: 6
    resources:
        runtime="6h",
    shell:
        """
        bcftools view -Oz {input.bcf} > {output.vcf}
        tabix {output.vcf}
        """


# rules for calling genotypes per individual (not used in final manuscript)


rule bcftools_mpileup_call_individual_chunk:
    """
    Call genotypes per sample in filtered sites regions. Filter for minimum
    depth, as well as heterozygous sites with poor allelic balance.
    """
    input:
        alignments="results/datasets/{dataset}/bams/{population}.{ref}{dp}.bam",
        bai="results/datasets/{dataset}/bams/{population}.{ref}{dp}.bam.bai",
        ref="results/ref/{ref}/{ref}.fa",
        index="results/ref/{ref}/{ref}.fa.fai",
        regions="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
        sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.bed",
    output:
        bcf=temp(
            "results/datasets/{dataset}/bcfs/chunks/{chunk}/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}_allbal{ablow}-{abhi}.bcf"
        ),
        idx=temp(
            "results/datasets/{dataset}/bcfs/chunks/{chunk}/{dataset}.{ref}_{population}{dp}_{sites}-filts.filtered_mindp{mindp}_allbal{ablow}-{abhi}.bcf.csi"
        ),
    conda:
        "../envs/bcftools121.yaml"
    threads: 2
    resources:
        runtime="1h",
    params:
        baseq=config["baseQ"],
        mapq=config["mapQ"],
    shell:
        """
        bcftools mpileup --threads {threads} -f {input.ref} -R {input.regions} \
            -Ou -T {input.sites} -B --min-MQ {params.mapq} \
            --min-BQ {params.baseq} -a "FORMAT/AD,FORMAT/DP,INFO/AD" \
            {input.alignments} | \
            bcftools call -m -f GQ,GP -Ou | \
            bcftools filter -g 5 -Ou | \
            bcftools view -V indels -Ou | \
            bcftools +setGT -Ou -- -t q -n . -i"FMT/DP<{wildcards.mindp}" | \
            bcftools +setGT -Ou -- -t q -n . \
                -i'GT="het" & (FMT/AD[:0]/FMT/DP < {wildcards.ablow} | FMT/AD[:0]/FMT/DP > {wildcards.abhi} | FMT/AD[:1]/FMT/DP < {wildcards.ablow} | FMT/AD[:1]/FMT/DP > {wildcards.abhi})' | \
            bcftools view -i'GT!="./."' -Ou | \
            bcftools +fill-tags -Ob -- -t all > {output.bcf}
        bcftools index -o {output.idx} {output.bcf}
        """


# function to get a list of BCFs for a grouping, here the whole dataset.
def bcfs(wildcards):
    dp = wildcards.dp
    samples = [
        sample
        for sample in angsd.samples.index.tolist()
        if sample not in config["calling_drop"]
    ]
    return {
        "bcfs": expand(
            "results/datasets/{{dataset}}/bcfs/chunks/{{chunk}}/{{dataset}}.{{ref}}_{population}{{dp}}_{{sites}}-filts.filtered_mindp{{mindp}}_allbal{{ablow}}-{{abhi}}.bcf",
            population=samples,
        ),
        "idxs": expand(
            "results/datasets/{{dataset}}/bcfs/chunks/{{chunk}}/{{dataset}}.{{ref}}_{population}{{dp}}_{{sites}}-filts.filtered_mindp{{mindp}}_allbal{{ablow}}-{{abhi}}.bcf.csi",
            population=samples,
        ),
    }


rule bcftools_merge_all:
    """
    Generate dataset-wide BCF from per sample BCFs. Par down to include
    monomorphic and biallelic sites and remove transition positions where DNA
    damage might be present. Sites are not yet filtered for missing data.
    """
    input:
        unpack(bcfs),
    output:
        bcf=temp(
            "results/datasets/{dataset}/bcfs/chunks/{chunk}/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-indcall_allbal{ablow}-{abhi}.notrans.bcf"
        ),
        idx=temp(
            "results/datasets/{dataset}/bcfs/chunks/{chunk}/{dataset}.{ref}_all{dp}_{sites}-filts.filtered_mindp{mindp}-allsites-indcall_allbal{ablow}-{abhi}.notrans.bcf.csi"
        ),
    conda:
        "../envs/bcftools121.yaml"
    threads: 2
    resources:
        runtime="2h",
    shell:
        """
        bcftools merge --force-samples -Ou {input.bcfs} | \
            bcftools +fill-tags -Ou -- -t all | \
            bcftools filter -g 5 -Ou | \
            bcftools view -M2 -V indels -Ov | \
            awk -F '\t' '!(($4 == "A" && $5 == "G") || ($4 == "G" && $5 == "A") || ($4 == "C" && $5 == "T") || ($4 == "T" && $5 == "C"))' | \
            bcftools view -Ob > {output.bcf}
        bcftools index -o {output.idx} {output.bcf}
        """
