//load file with paths to necessary software (meerkat paths)
load '/group/bioi1/jovanam/scripts/bpipe-pipelines/software_paths.groovy'
//load bpipe stages required for this pipeline
load '/group/bioi1/jovanam/scripts/bpipe-pipelines/rnaseq_stages.groovy'
load '/group/bioi1/jovanam/scripts/bpipe-pipelines/common_stages.groovy'

//path to adapter file for adapter type used in this experiment
ADAPTERS = "/usr/local/installed/trimmomatic/0.35/adapters/NexteraPE-PE.fa"
//trimmomatic parameters
LEADING = 25
TRAILING = 25
MINLEN = 30
ILLUMINACLIP = "2:30:10"

GEN_DIR = "\$GENOMES/hg38/star_2_5/" // directory of the genome index to be used for mapping
GEN_FASTA = "\$GENOMES/hg38/fasta/hg38.fa" // genome fasta file
GTF = "\$GENOMES/hg38/star/hg38_GENCODEV20_Comp.gtf" // GTF file (used by STAR for junction mapping)
SAF = "\$GENOMES/hg38/saf/hg38_GENCODEV20_Comp.saf" // SAF file (used by featureCounts for read counting)
SJBOHANG = 75 //STAR: length of genomic sequence around annotated junctioin; max(ReadLength)-1; default: 100

//number of treads per sample for multithreaded tools
NTHREADS = 4 // no. threads to use for mapping, then counting

MULTIQCDIR = "."

//paired-end rna-seq pipeline, starting with data from SRA
run {
    "%.fastq.gz" * [ fastqc ] +
    "%_R*.fastq.gz" * [ trim_PE ] + star_map_1pass_PE +
    "%_R*.fastq.trim.gz" * [ star_map_2pass_PE ] + count_reads_RNA + 
    "%.bam" * [ sort_bam + index_bam ] + multiqc 
}

