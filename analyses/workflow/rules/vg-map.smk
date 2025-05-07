"""
This Snakefile contains all the rules related to mapping samples to the variant
graph made from the species reference and modern sample vcfs. This graph has to
first be made using a separate snakemake run. After it is made, a Snakemake run
with both modern and historical samples will be run that maps them both to the
new variant graph instead of the linear reference. All the rules for that
process are here. Several are imported from PopGLen to reuse their setup, while
changing outputs so that both bwa and vg mapped bams can live together in the
same working directory for testing and comparison.
"""


use rule fastp_mergedout from angsd as vg_fastp_mergedout with:
    # Trim raw reads with fastp, merging overlapping reads (for historical
    # samples). Runs the same as PopGLen, but makes a separate storage location
    # to avoid conflicts with files from default PopGLen runs in the same
    # folder.
    output:
        trimmed=temp(
            expand(
                "results/preprocessing/fastp/vg/{{sample}}_{{unit}}_{{lib}}.{read}.uncollapsed.fastq.gz",
                read=["R1", "R2"],
            )
        ),
        merged=temp(
            "results/preprocessing/fastp/vg/{sample}_{unit}_{lib}.merged.fastq.gz"
        ),
        html=report(
            "results/preprocessing/qc/fastp/vg/{sample}_{unit}_{lib}.merged.html",
            category="00 Quality Control",
            subcategory="1 Trimming Reports (for vg)",
            labels={
                "Sample": "{sample}",
                "Unit": "{unit}",
                "Lib": "{lib}",
                "Type": "fastp Report",
            },
        ),
        json="results/preprocessing/qc/fastp/vg/{sample}_{unit}_{lib}.merged.json",
    log:
        "logs/preprocessing/fastp/vg/{sample}_{unit}_{lib}.merged.log",
    benchmark:
        "benchmarks/preprocessing/fastp/vg/{sample}_{unit}_{lib}.merged.log"
    threads: lambda w, attempt: attempt * 12


use rule fastp_pairedout from angsd as vg_fastp_pairedout with:
    # Trim raw reads with fastp, merging overlapping reads (for historical
    # samples). Runs the same as PopGLen, but makes a separate storage location
    # to avoid conflicts with files from default PopGLen runs in the same
    # folder.
    output:
        trimmed=temp(
            expand(
                "results/preprocessing/fastp/vg/{{sample}}_{{unit}}_{{lib}}.{read}.paired.fastq.gz",
                read=["R1", "R2"],
            )
        ),
        html=report(
            "results/preprocessing/qc/fastp/vg/{sample}_{unit}_{lib}.paired.html",
            category="00 Quality Control",
            subcategory="1 Trimming Reports (for vg)",
            labels={
                "Sample": "{sample}",
                "Unit": "{unit}",
                "Lib": "{lib}",
                "Type": "fastp Report",
            },
        ),
        json="results/preprocessing/qc/fastp/vg/{sample}_{unit}_{lib}.paired.json",
    log:
        "logs/preprocessing/fastp/vg/{sample}_{unit}_{lib}.paired.log",
    benchmark:
        "benchmarks/preprocessing/fastp/vg/{sample}_{unit}_{lib}.paired.log"
    threads: lambda w, attempt: attempt * 12


rule vg_map_merged:
    """
    Maps merged fastq files to the variant graph. Adds a read group, the MD tag
    with the reference states.
    """
    input:
        fq="results/preprocessing/fastp/vg/{sample}_{unit}_{lib}.merged.fastq.gz",
        fa="results/ref/{ref}/{ref}.fa",
        fai="results/ref/{ref}/{ref}.fa.fai",
        vg="results/ref/{ref}/{ref}.modvars.vg",
        xg="results/ref/{ref}/{ref}.modvars.xg",
        gcsa="results/ref/{ref}/{ref}.modvars.gcsa",
    output:
        bam=temp("results/mapping/mapped/{sample}_{unit}_{lib}.{ref}.vg.merged.bam"),
        flagstat="results/mapping/mapped/{sample}_{unit}_{lib}.{ref}.vg.merged.bam.flagstat",
    conda:
        "../envs/vg.yaml"
    params:
        rg=angsd.get_read_group,
    threads: 64
    resources:
        runtime="1d",
    shell:
        """
        vg map -t {threads} -w 300 -k 15 --log-time \
            -x {input.xg} \
            -g {input.gcsa} \
            -f {input.fq} \
            --surject-to bam | \
            samtools addreplacerg -u -r {params.rg} - | \
            samtools calmd -u - {input.fa} | \
            samtools sort -T {resources.tmpdir} -o {output.bam}
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        """


rule vg_map_paired:
    """
    Maps paired fastq files to the variant graph. Adds a read group, the MD tag
    with the reference states.
    """
    input:
        fq1="results/preprocessing/fastp/vg/{sample}_{unit}_{lib}.R1.paired.fastq.gz",
        fq2="results/preprocessing/fastp/vg/{sample}_{unit}_{lib}.R2.paired.fastq.gz",
        fa="results/ref/{ref}/{ref}.fa",
        fai="results/ref/{ref}/{ref}.fa.fai",
        vg="results/ref/{ref}/{ref}.modvars.vg",
        xg="results/ref/{ref}/{ref}.modvars.xg",
        gcsa="results/ref/{ref}/{ref}.modvars.gcsa",
    output:
        bam=temp("results/mapping/mapped/{sample}_{unit}_{lib}.{ref}.vg.paired.bam"),
        flagstat="results/mapping/mapped/{sample}_{unit}_{lib}.{ref}.vg.paired.bam.flagstat",
    conda:
        "../envs/vg.yaml"
    params:
        rg=angsd.get_read_group,
    threads: 64
    resources:
        runtime="1d",
    shell:
        """
        vg map -t {threads} --log-time \
            -x {input.xg} \
            -g {input.gcsa} \
            -f {input.fq1} -f {input.fq2} \
            --surject-to bam | \
            samtools addreplacerg -u -r {params.rg} - | \
            samtools calmd -u - {input.fa} | \
            samtools sort -T {resources.tmpdir} -o {output.bam}
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        """


rule picard_reorder:
    """
    Uses picard to sort everything in the mapped BAMs in the same order as the
    reference fasta, which vg doesn't do automatically.
    """
    input:
        bam="results/mapping/mapped/{sample}_{unit}_{lib}.{ref}.vg.{pairing}.bam",
        fa="results/ref/{ref}/{ref}.fa",
        fai="results/ref/{ref}/{ref}.fa.fai",
    output:
        bam=temp(
            "results/mapping/mapped/{sample}_{unit}_{lib}.{ref}.vg.{pairing}.reordered.bam"
        ),
    conda:
        "../envs/vg.yaml"
    threads: 32
    resources:
        runtime="6h",
    shell:
        """
        picard ReorderSam -Xmx{resources.mem_mb}m --INPUT {input.bam} \
                --SEQUENCE_DICTIONARY {input.fa} --OUTPUT {output.bam}
        """


def get_vg_unit_bams(wildcards):
    reads = angsd.units.loc[angsd.units["sample"] == wildcards.sample]
    combos = reads[["sample", "unit", "lib"]].agg("_".join, axis=1)
    return expand(
        "results/mapping/mapped/{combo}.{{ref}}.vg.{{pairing}}.reordered.bam",
        combo=combos,
    )


use rule samtools_merge_paired_units from angsd as samtools_merge_units_vg with:
    # Merges sequencing runs per sample (in this study, no samples had multiple
    # libraries) after they've been mapped with vg.
    input:
        get_vg_unit_bams,
    output:
        bam=temp("results/mapping/mapped/{sample}.{ref}.vg.{pairing}.bam"),
    log:
        "logs/mapping/samtools/merge/{sample}.{ref}.vg.{pairing}.log",
    benchmark:
        "benchmarks/mapping/samtools/{sample}.{ref}.vg.{pairing}.log"


use rule mark_duplicates from angsd as mark_duplicates_vg with:
    # Remove duplicate reads from paired end bam files mapped with vg
    input:
        bams="results/mapping/mapped/{sample}.{ref}.vg.paired.bam",
        flagstat="results/mapping/mapped/{sample}.{ref}.vg.paired.flagstat",
    output:
        bam=temp("results/mapping/dedup/{sample}.{ref}.vg.paired.rmdup.bam"),
        metrics="results/mapping/qc/mark_duplicates/{sample}.{ref}.vg.paired.picard.metrics",
    log:
        "logs/mapping/picard/dedup/{sample}.{ref}.vg.paired.log",
    benchmark:
        "benchmarks/mapping/picard/dedup/{sample}.{ref}.vg.paired.log"
    threads: 10


rule dedup_merged_vg:
    # Remove duplicates from collapsed read bam files mapped with vg
    input:
        bam="results/mapping/mapped/{sample}.{ref}.vg.merged.bam",
        flagstat="results/mapping/mapped/{sample}.{ref}.vg.merged.flagstat",
    output:
        json="results/mapping/qc/dedup/{sample}.{ref}.vg.merged.dedup.json",
        hist="results/mapping/qc/dedup/{sample}.{ref}.vg.merged.dedup.hist",
        log="results/mapping/qc/dedup/{sample}.{ref}.vg.merged.dedup.log",
        bam=temp("results/mapping/dedup/{sample}.{ref}.vg.merged_rmdup.bam"),
        bamfin="results/mapping/bams/{sample}.{ref}.vg.merged.rmdup.bam",
        baifin="results/mapping/bams/{sample}.{ref}.vg.merged.rmdup.bam.bai",
    conda:
        "../envs/dedup.yaml"
    shadow:
        "minimal"
    threads: lambda wildcards, attempt: attempt * 8
    params:
        outdir=lambda w, output: os.path.dirname(output.bam),
    resources:
        runtime=lambda wildcards, attempt: attempt * 1440,
    shell:
        """
        dedup -i {input.bam} -m -u -o {params.outdir}
        samtools sort -T {resources.tmpdir} -o {output.bamfin} {output.bam}
        samtools index {output.bamfin}
        mv {params.outdir}/{wildcards.sample}.{wildcards.ref}.vg.merged.dedup.json \
            {output.json}
        mv {params.outdir}/{wildcards.sample}.{wildcards.ref}.vg.merged.hist \
            {output.hist}
        mv {params.outdir}/{wildcards.sample}.{wildcards.ref}.vg.merged.log \
            {output.log}
        """


use rule bam_clipoverlap from angsd as bam_clipoverlap_vg with:
    # Clip overlapping reads in paired end bam files mapped with vg
    input:
        bam="results/mapping/dedup/{sample}.{ref}.vg.paired.rmdup.bam",
        ref="results/ref/{ref}/{ref}.fa",
    output:
        bam=temp("results/mapping/dedup/{sample}.{ref}.vg.paired.rmdup.clip.bam"),
        log="results/mapping/qc/bamutil_clipoverlap/{sample}.{ref}.vg.paired.rmdup.clipoverlap.stats",
    log:
        "logs/mapping/bamutil/clipoverlap/{sample}.{ref}.vg.paired.rmdup.log",
    benchmark:
        "benchmarks/mapping/bamutil/clipoverlap/{sample}.{ref}.vg.paired.rmdup.log"
    threads: 6


rule finalize_paired_vg_bam:
    """
    Move the clipped paired end bam file from the vg mapping rules to its final
    location and index it.
    """
    input:
        "results/mapping/dedup/{sample}.{ref}.vg.paired.rmdup.clip.bam",
    output:
        bam="results/mapping/bams/{sample}.{ref}.vg.paired.rmdup.clip.bam",
        bai="results/mapping/bams/{sample}.{ref}.vg.paired.rmdup.clip.bam.bai",
    container:
        angsd.samtools_container
    shell:
        """
        cp {input} {output.bam}
        samtools index {output.bam}
        """
