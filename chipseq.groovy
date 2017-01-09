//load file with paths to necessary software (meerkati paths)
load '/group/bioi1/jovanam/scripts/bpipe-pipelines/software_paths.groovy'
//load bpipe stages required for this pipeline
load '/group/bioi1/jovanam/scripts/bpipe-pipelines/dnaseq_stages.groovy'
load '/group/bioi1/jovanam/scripts/bpipe-pipelines/common_stages.groovy'

//path to adapter file for adapter type used in this experiment
ADAPTERS = "/usr/local/installed/trimmomatic/0.35/adapters/TruSeq2-SE.fa"
//trimmomatic parameters
LEADING = 25
TRAILING = 25
MINLEN = 30
ILLUMINACLIP = "2:30:10"

//path to bowtie2 index
INDEX = "\$GENOMES/hg19/bowtie2-index/hg19"
//number of treads per sample for multithreaded tools
NTHREADS = 4

//single-end chip-seq pipeline, starting with data from SRA
run {
  "%.sra" * [ sra_to_fastq_SE + fastqc + trim_SE + bowtie2_map_SE ] +
  "%.bam" * [ sort_bam + index_bam + dedup ] + multiqc
}

