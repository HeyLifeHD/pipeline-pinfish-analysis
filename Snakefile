## General pipeline parameters:
# ABSOLUTE path to directory holding the working directory:
workdir_top: "/bigdisk/Nanopore/"

# Name of the pipeline:
pipeline: "190220_pinfish_analysis"

# Repository URL:
repo: "https://github.com/HeyLifeHD/pipeline-pinfish-analysis.git"

## Pipeline-specific parameters:
# Number of threads:
threads: 5

# Input genome:
genome_fasta: "/home/epicwl/genomes/hg19/male.hg19.fa"

# CDNA or direct RNA reads in fastq format:
reads_fastq: "/bigdisk/Nanopore/raw_data/DAC_SB939.fastq"

# Decide if pychopper should be used and give path to adapter sequences
pychopper_analysis: true
adapters_fastq: "/home/epicwl/pychopper/data/cdna_barcodes.fas"

#If no pychopper is used decide which phred score and length cut off should be used for fastq filtering. In addition specify hard headcrop if needed
min_phred: "7"
min_length: "500"
headcrop: "75"

# Extra option passed to minimap2 when generating index:
#kmer length:
minimap_index_opts: "-k14"

## Extra options passed to minimap2 when mapping:
# Enable this for stranded data (this is the case if pychopper was run):
minimap2_opts: "-uf"

# Enable this for SIRV data:
#minimap2_opts: "--splice-flank=no"

# Minmum mapping quality:
minimum_mapping_quality: 10

# Options passed to spliced_bam2gff:
# Enable this for stranded data (this is the case if pychopper was run):
spliced_bam2gff_opts: "-s"

# -c parameter:
minimum_cluster_size: 5

# -p parameter:
minimum_isoform_percent: 1.0

# -d parameter:
exon_boundary_tolerance: 15

# -e parameter:
terminal_exon_boundary_tolerance: 45

# Extra options passed to minimap2 when mapping polished reads:
# Enable this for stranded data (this is the case if pychopper was run):
minimap2_opts_polished: "-uf"
# Enable this for SIRV data:
#minimap2_opts_polished: "--splice-flank=no"

# Options passed to spliced_bam2gff when converting alignments of polished reads:
# Enable this for stranded data (this is the case if pychopper was run):
spliced_bam2gff_opts_pol: "-s"

# Options passed to collapse_partials when collapsing fragmentation artifacts
# in clustered and polished transcripts.

# Internal exon boundary tolerance:
collapse_internal_tol: 7

# Five prime boundary tolerance:
collapse_five_tol: 7500

# Three prime boundary tolerance:
collapse_three_tol: 30