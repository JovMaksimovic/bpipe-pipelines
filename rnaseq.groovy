load 'rnaseq_stages.groovy'
load 'common_stages.groovy'

ADAPTERS = "/usr/local/Trimmomatic-0.33/adapters/TruSeq2-SE.fa" // file of adapter sequence to be used by Trimmomatic

GEN_DIR = "/mnt/storage/shared/genomes/hg19/star" // directory of the genome to be used for mapping
GEN_FASTA = "/mnt/storage/shared/genomes/hg19/fasta/hg19.fa" // genome fasta file
GTF = "/mnt/storage/shared/genomes/hg19/star/hg19_gencodeV19_comp.gtf" // GTF file (used by STAR for junction mapping)

SAF = "/mnt/storage/shared/genomes/hg19/saf/hg19_GENCODEV19_Comp.saf" // SAF file (used by featureCounts for read counting)

MAPPING_DIR = "mapped" // name of directory where mapped reads will go
NTHREADS = 20 // no. threads to use for mapping, then counting

run {
    // ensure that filenames suffixes match *actual* filenames
    "%.fastq.gz" * [ fastqc + trim] + star_map_1pass_SE + star_gen2pass +
    "%_L00*_R1_001.fastq.trim.gz" * [ star_map_2pass_SE + sort_bam + index_bam ] + count_reads_RNA
}

