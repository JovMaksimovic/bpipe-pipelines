sra_to_fastq_PE = {
  // Convert SRA files to fastq format, paired-end 
  output.dir = "fastq"
  transform(".sra") to ("_1.fastq.gz","_2.fastq.gz"){
    exec """
      $FASTQ_DUMP --split-3 --gzip -O $output.dir $input.sra;
    ""","small"
  }
}

sra_to_fastq_SE = {
  // Convert SRA files to fastq format, single-end
  output.dir = "fastq"
  transform(".sra") to (".fastq.gz"){
    exec """
      $FASTQ_DUMP --split-3 --gzip -O $output.dir $input.sra;
    ""","small"
  }
}

concat_fastq_SE = {
    // Concatenate fastq files from the same sample
    doc "Concatenate fastq files from the same run"
    output.dir = "catfastq"

    produce(output.gz.prefix.prefix.prefix + ".cat.fastq.gz"){
        exec """
                cat $inputs.gz > $output.gz
        ""","tiny"
    }
}

fastqc = {
    // fastqc quality control
    doc "Quality control using FASTQC"
    output.dir = "fastqc"

    transform(".fastq.gz") to ("_fastqc.zip") {
        exec """
		$FASTQC -o $output.dir $input.gz
	""","tiny"
    }
    println inputs;
    forward inputs
}

multiqc = {
   // summarise statistics from all tools using multiqc
   doc "Pipeline summary report by MultiQC"

   exec """
	multiqc . --ignore .bpipe
   ""","multiqc"
}

trim_SE = {
    // trim single-end reads using trimmomatic
    doc "Trim poor quility bases and/or adapter sequences from reads using Trimmomatic"
    output.dir = "trimmed"

    filter("trim"){
        exec """ 
            $TRIMOMMATIC SE -threads $NTHREADS $input.gz $output.gz ILLUMINACLIP:$ADAPTERS:$ILLUMINACLIP 
            TRAILING:$TRAILING LEADING:$LEADING MINLEN:$MINLEN
        ""","trimmomatic"
    }
}

trim_PE = {
   // trim paired-end reads using trimmomatic
   doc "Trim poor quility bases and/or adapter sequences from reads using Trimmomatic"
   output.dir="trimmed"

   filter("trim","trim") {
        exec """
                $TRIMOMMATIC PE -threads $NTHREADS
                $input1.gz $input2.gz
                $output1.gz ${output1.prefix}.unpaired.gz
                $output2.gz ${output2.prefix}.unpaired.gz
		ILLUMINACLIP:$ADAPTERS:$ILLUMINACLIP
                LEADING:$LEADING TRAILING:$TRAILING MINLEN:$MINLEN
        ""","trimmomatic"
   }
}

sort_bam = {
    // sort bam files
    doc "Sort bam files"
    output.dir = "sorted"

    filter("srt") {
        exec """
		$SAMTOOLS sort -m 1G -@ $NTHREADS -o $output $input
	""","srtindex"
    }
}

index_bam = {
    // index bam files
    doc "Index bam files"
    output.dir="sorted"

    transform("bam") to ("bam.bai") {
        exec """
		$SAMTOOLS index $input.bam
	""","srtindex"
    }
    forward input
}

count_mapped_SE = {
    // Count total number of single-end mapped reads in bam file (no multi-map)
    doc "Count total number of single-end mapped reads in bam file (no multi-map)"
    output.dir = "stats"

    transform("txt"){
        exec """
                $SAMTOOLS view -F 0x904 -c $input > $output
        ""","tiny"
    }
}

count_mapped_PE = {
    // Count total number of mapped paired-end reads in bam file (no multi-map)
    doc "Count total number of mapped paired-end reads in bam file (no multi-map)"
    output.dir = "stats"

    transform("txt"){
        exec """
		$SAMTOOLS view -F 0x4 $input.bam | cut -f 1 | sort | uniq | wc -l > $output
        ""","tiny"
    }
}

coverage_by_region = {
    // Get read coverage for each base across a particular region(s)
    doc "Get read coverage for each base across a particular region(s)"
    output.dir = "stats"

    transform("bed"){
        exec """
                $BEDTOOLS coverage -d -split -abam $input.bam -b $REGIONS_BED > $output.bed
        ""","regioncov"
    }
}

