//load file with paths to necessary software (meerkat paths)
load '/group/bioi1/jovanam/scripts/bpipe-pipelines/software_paths.groovy'
//load bpipe stages required for this pipeline
load '/group/bioi1/jovanam/scripts/bpipe-pipelines/dnaseq_stages.groovy'
load '/group/bioi1/jovanam/scripts/bpipe-pipelines/common_stages.groovy'

//path to adapter file for adapter type used in this experiment
ADAPTERS = "/usr/local/installed/trimmomatic/0.35/adapters/NexteraPE-PE.fa"
//trimmomatic parameters
LEADING = 25
TRAILING = 25
MINLEN = 30
ILLUMINACLIP = "2:30:10" 

//path to bowtie2 index
INDEX = "\$GENOMES/hg19/bowtie2-index/hg19"
//number of treads per sample for multithreaded tools
NTHREADS = 4

//paired-end atac-seq pipeline, starting with data from SRA
run {
  "%.sra" * [ sra_to_fastq_PE + fastqc ] +
  "%_*.fastq.gz" * [ trim_PE + bowtie2_map_PE ] u+
  "%.bam" * [ sort_bam + index_bam + dedup ] + multiqc
}

