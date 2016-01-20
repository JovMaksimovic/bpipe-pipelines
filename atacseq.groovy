load 'dnaseq_stages.groovy'
load 'common_stages.groovy'

ADAPTERS_FASTA = "/usr/local/Trimmomatic-0.33/adapters/TruSeq2-PE.fa"
INDEX = "/mnt/storage/shared/genomes/hg19/bowtie2-index/hg19"
SAF = "/mnt/BLOO1-NGS/ATAC_seq/data/refGeneWithProm.saf"
NTHREADS = 10

run {
  "%.fastq.gz" * [fastqc] +
  "%_*.fastq.gz" * [trim_PE + bowtie2_map_PE] +
  "%.bam" * [sort_bam + index_bam + dedup] + count_reads_DNA
}

