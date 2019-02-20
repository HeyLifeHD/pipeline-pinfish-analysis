
import os
from os import path

#Directory of Snakefile for pinfish paths
SNAKEDIR = path.dirname(workflow.snakefile)

#Specify configfile
configfile: path.join(SNAKEDIR, "config.yaml")
workdir: path.join(config["workdir_top"], config["pipeline"])

#Additional Rules
include: "snakelib/utils.snake"


rule build_minimap_index: ## build minimap2 index
    input:
        genome = config["genome_fasta"]
    output:
        index = "index/genome_index.mmi"
    params:
        opts = config["minimap_index_opts"]
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
        minimap2 -t {threads} {params.opts} -I 1000G -d {output.index} {input.genome}
    """

rule raw_qc: # do qc of raw fastq
    input:
        fastq = config["reads_fastq"]
    output:
        txt = "QC/raw/NanoStats.txt"
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
        NanoPlot -t {threads} --fastq {input.fastq} -o QC/raw/
        touch {output.txt}
    """

rule read_directionality:
    input:
        fastq = config["reads_fastq"],
        adapters = config["adapters_fastq"]
    output:
        report = "pychopper/report.pdf",
        unclassified = "pychopper/unclassified.fq",
        classified = "pychopper/classified.fq"
    conda: "env.yml"
    shell:"""
        cdna_classifier.py -b {input.adapters] -r {output.report} \
        -u {output.unclassified} {input.fastq} {output.classified}
    """

rule pychopper_qc: # do qc of raw fastq
    input:
        fastq = "pychopper/classified.fq"
    output:
        txt = "QC/pychopper/NanoStats.txt"
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
        NanoPlot -t {threads} --fastq {input.fastq} -o QC/pychopper/
        touch {output.txt}
    """

rule filter_rawfastq:
    input:
        index = "index/genome_index.mmi",
        fastq = config["reads_fastq"]
    output:
        fastq = "filtered/filtered.fastq"
    params:
        qual = config["min_phred"],
        length = config["min_length"],
        head = config["headcrop"],
    threads: config["threads"]
    conda: "env.yml"
    shell:"""
        cat {input.fastq} \
        | NanoFilt -q {params.qual} -l {params.length} --headcrop {params.head} \
        > {output.fastq}
    """       

rule filter_qc: # do qc of raw fastq
    input:
        fastq = "filtered/filtered.fastq"
    output:
        txt = "QC/filter/NanoStats.txt"
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
        NanoPlot -t {threads} --fastq {input.fastq} -o QC/raw/
        touch {output.txt}
    """

rule map_reads: ## map reads using minimap2
    input:
       index = "index/genome_index.mmi",
       fastq = if config["pychopper_analysis"] "pychopper/classified.fq" else "filtered/filtered.fastq"
    output:
       bam = "alignments/reads_aln_sorted.bam"
    params:
        opts = config["minimap2_opts"],
        min_mq = config["minimum_mapping_quality"]
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
        minimap2 -t {threads} -ax splice {params.opts} {input.index} {input.fastq}\
        | samtools view -q {params.min_mq} -F 2304 -Sb | samtools sort -@ {threads} - -o {output.bam};
        samtools index {output.bam}
    """

rule bam_qc: # do qc of raw fastq
    input:
        bam = "alignments/reads_aln_sorted.bam"
    output:
        txt = "QC/bam/NanoStats.txt"
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
        NanoPlot -t {threads} --bam {input.bam} -o QC/bam/
        touch {output.txt}
    """

rule convert_bam: ## convert BAM to GFF
    input:
        bam = "alignments/reads_aln_sorted.bam"
    output:
        raw_gff = "results/raw_transcripts.gff"
    conda: "env.yml"
    params:
        opts = config["spliced_bam2gff_opts"]
    threads: config["threads"]
    shell:"""
    {SNAKEDIR}/pinfish/spliced_bam2gff/spliced_bam2gff {params.opts} -t {threads} -M {input.bam} > {output.raw_gff}
    """

rule cluster_gff: ## cluster transcripts in GFF
    input:
        raw_gff = "results/raw_transcripts.gff"
    output:
        cls_gff = "results/clustered_transcripts.gff",
        cls_tab = "results/cluster_memberships.tsv",
    params:
        c = config["minimum_cluster_size"],
        d = config["exon_boundary_tolerance"],
        e = config["terminal_exon_boundary_tolerance"],
        min_iso_frac = config["minimum_isoform_percent"],
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
    {SNAKEDIR}/pinfish/cluster_gff/cluster_gff -p {params.min_iso_frac} -t {threads} -c {params.c} -d {params.d} -e {params.e} -a {output.cls_tab} {input.raw_gff} > {output.cls_gff}
    """

rule collapse_clustered: ## collapse clustered read artifacts
    input:
        cls_gff = "results/clustered_transcripts.gff"
    output:
        cls_gff_col = "results/clustered_transcripts_collapsed.gff"
    params:
        d = config["collapse_internal_tol"],
        e = config["collapse_three_tol"],
        f = config["collapse_five_tol"],
    conda: "env.yml"
    shell:"""
    {SNAKEDIR}/pinfish/collapse_partials/collapse_partials -d {params.d} -e {params.e} -f {params.f} {input.cls_gff} > {output.cls_gff_col}
    """

rule polish_clusters: ## polish read clusters
    input:
        cls_gff = "results/clustered_transcripts.gff",
        cls_tab = "results/cluster_memberships.tsv",
        bam = "alignments/reads_aln_sorted.bam"
    output:
        pol_trs = "results/polished_transcripts.fas",
    params:
        c = config["minimum_cluster_size"],
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
    {SNAKEDIR}/pinfish/polish_clusters/polish_clusters -t {threads} -a {input.cls_tab} -c {params.c} -o {output.pol_trs} {input.bam}
    """

rule map_polished: ## map polished transcripts to genome
    input:
       index = "index/genome_index.mmi",
       fasta = "results/polished_transcripts.fas"
    output:
       pol_bam = "alignments/polished_reads_aln_sorted.bam"
    params:
        extra = config["minimap2_opts_polished"]
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} {params.extra} -ax splice {input.index} {input.fasta}\
    | samtools view -Sb -F 2304 | samtools sort -@ {threads} - -o {output.pol_bam};
    samtools index {output.pol_bam}
    """

rule convert_polished: ## convert BAM of polished transcripts to GFF
    input:
        bam = "alignments/polished_reads_aln_sorted.bam"
    output:
        pol_gff = "results/polished_transcripts.gff"
    params:
        extra = config["spliced_bam2gff_opts_pol"]
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
    {SNAKEDIR}/pinfish/spliced_bam2gff/spliced_bam2gff {params.extra} -t {threads} -M {input.bam} > {output.pol_gff}
    """

rule collapse_polished: ## collapse polished read artifacts
    input:
        pol_gff = "results/polished_transcripts.gff"
    output:
        pol_gff_col = "results/polished_transcripts_collapsed.gff"
    params:
        d = config["collapse_internal_tol"],
        e = config["collapse_three_tol"],
        f = config["collapse_five_tol"]
    conda: "env.yml"
    shell:"""
    {SNAKEDIR}/pinfish/collapse_partials/collapse_partials -d {params.d} -e {params.e} -f {params.f} {input.pol_gff} > {output.pol_gff_col}
    """

rule gen_corr_trs: ## Generate corrected transcriptome.
    input:
        genome = config["genome_fasta"],
        gff = "results/polished_transcripts_collapsed.gff"
    output:
        fasta = "results/corrected_transcriptome_polished_collapsed.fas"
    conda: "env.yml"
    shell:"""
    gffread -g {input.genome} -w {output.fasta} {input.gff}
    """

rule compare_gff: ## compare gff with reference
    input:
        reference = config["reference_transcriptome"],
        gff = "results/polished_transcripts_collapsed.gff"
    output:
        stats = "results/compare/compare.stats"
    conda: "env.yml"
    shell:"""
    gffread {input.reference} -T results/polished_transcripts_collapsed.gtf
    gffcompare -R -r results/polished_transcripts_collapsed.gtf -o results/compare/compare {input.gff}
    touch {output.stats}
    """

rule all: ## run the whole pipeline
    input:
        index = rules.build_minimap_index.output.index,
        raw_qc = rules.raw_qc.output.txt,
        if config["pychopper_analysis"] "QC/pychopper/NanoStats.txt" else "QC/filter/NanoStats.txt",
        aligned_reads = rules.map_reads.output.bam,
        bam_qc = rules.bam_qc.output.txt,
        raw_gff = rules.convert_bam.output.raw_gff,
        cls_gff = rules.cluster_gff.output.cls_gff,
        cls_gff_col = rules.collapse_clustered.output.cls_gff_col,
        pol_trs = rules.polish_clusters.output.pol_trs,
        pol_bam = rules.map_polished.output.pol_bam,
        pol_gff = rules.convert_polished.output.pol_gff,
        pol_gff_col = rules.collapse_polished.output.pol_gff_col,
        corr_trs = rules.gen_corr_trs.output.fasta,
        gff_compare = rules.compare_gff.ouput.stats,
