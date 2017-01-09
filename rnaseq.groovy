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

GEN_DIR = "\$GENOMES/hg19/star-new" // directory of the genome index to be used for mapping
GEN_FASTA = "\$GENOMES/hg19/fasta/hg19.fa" // genome fasta file
GTF = "\$GENOMES/hg19/star/hg19_gencodeV19_comp.gtf" // GTF file (used by STAR for junction mapping)
SAF = "\$GENOMES/hg19/saf/hg19_GENCODEV19_Comp.saf" // SAF file (used by featureCounts for read counting)
SJBOHANG = 100 //STAR: length of the genomic sequence around the annotated junctioin; max(ReadLength)-1; default: 100

//number of treads per sample for multithreaded tools
NTHREADS = 4 // no. threads to use for mapping, then counting

//paired-end rna-seq pipeline, starting with data from SRA
run {
    "%.sra" * [ sra_to_fastq_PE + fastqc ] +
    "%_*.fastq.gz" * [ trim_PE ] + star_map_1pass_PE +
    "%_*.fastq.trim.gz" * [ star_map_2pass_PE + sort_bam + index_bam ] + count_reads_RNA + multiqc 
}

