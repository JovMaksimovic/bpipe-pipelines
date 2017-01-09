bowtie2_map_SE = {
  //Map single-end reads using bowtie2
  doc "Map single-end reads using bowtie2"
  String lanes = "${inputs}";

  output.dir = "mapped"
  transform("bam"){
    exec """
       $BOWTIE2 -x $INDEX -p $NTHREADS --local --no-unal -U $input.gz | samtools view -S -b - > $output.bam
    ""","bowtie2"
  }
}

bowtie2_map_PE = {
    //String files = inputs;
    //def allFiles = files.tokenize(' ');
    //R1 = allFiles.findAll{ item -> item.contains('R1') }.join(',')
    //R2 = allFiles.findAll{ item -> item.contains('R2') }.join(',')

    //Map paired-end reads using bowtie2
    doc "Map paired-end reads using bowtie2"
    output.dir = "mapped"

    transform("bam"){
        exec """
            $BOWTIE2 -x $INDEX -p $NTHREADS --local --no-unal -1 $input1.gz -2 $input2.gz | samtools view -S -b - > $output.bam
        ""","bowtie2"
    }
}

bowtie_map_SE = {
    //Map single-end reads using bowtie
    doc "Map single-end reads using bowtie"

    output.dir = "mapped"
    transform("bam"){
        exec """
            $BOWTIE --sam -p $NTHREADS -q $INDEX $input.gz | samtools view -S -b - > $output.bam
        ""","bowtie"
    }
}

dedup = {
    //Remove PCR duplicates from reads
    doc "Remove PCR duplicates from reads"
    output.dir="deduped"

        exec """
            java -jar -Xmx4g $MARK_DUPLICATES
                INPUT=$input.bam
                REMOVE_DUPLICATES=true
                VALIDATION_STRINGENCY=LENIENT
                AS=true
                METRICS_FILE=$output.metrics
                CREATE_INDEX=true
                OUTPUT=$output.bam
        ""","small"
}

count_reads_DNA = {
    //Count reads across features using DNA data
    doc "Count reads across features using DNA"
    output.dir = "counts"

    produce("counts.txt") {
        exec """
            $FEATURECOUNTS -O -g GeneID -T $NTHREADS -F SAF -a $SAF -o $output $inputs.bam
        ""","count"
    }
}

