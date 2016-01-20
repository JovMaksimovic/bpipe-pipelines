load 'dnaseq_stages.groovy'
load 'common_stages.groovy'

ADAPTERS_FASTA = "/usr/local/Trimmomatic-0.33/adapters/TruSeq2-SE.fa"
INDEX = "/mnt/storage/shared/genomes/mm10/bowtie2-index/mm10"
SAF = "/mnt/CELL1-NGS/ChIPseq_RNAseq_Betty_Kao/ChIPseqHistoneMeth_Analysis/mm10GenesUnique.saf"
NTHREADS = 10

run {
    "%.fastq.gz" * [ fastqc + trim] +
    ~"(.*)_ML150.*_ML150371_15090[89]_HG.*BGXX_NS500295_L00[1-4]_R1.fastq.trim.gz" * [ bowtie2_map + sort_bam + index_bam] + 
    count_reads
}

