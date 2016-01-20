bowtie2_map = {
  String lanes = "${inputs}";

  lanes = lanes.replaceAll(' ',',');

  output.dir = "mapped"
  transform("bam"){
    exec """
       bowtie2 -x $INDEX -p $NTHREADS --no-unal -U $lanes | $SAMTOOLS view -S -b - > $output.bam
    """
  }
}

bowtie2_map_PE = {

    String files = inputs;
    def allFiles = files.tokenize(' ');
    R1 = allFiles.findAll{ item -> item.contains('R1') }.join(',')
    R2 = allFiles.findAll{ item -> item.contains('R2') }.join(',')

    output.dir = "mapped"
    transform("bam"){
        exec """
            bowtie2 -x $INDEX -p $NTHREADS --no-unal -1 ${R1} -2 ${R2} | samtools view -S -b - > $output.bam
        """
    }
}

dedup = {
    doc "Remove PCR duplicates from reads"
    output.dir="deduped"

    transform("dedup.bam",".metrics"){
        exec """
            java -jar -Xmx4g $PICARD_HOME/lib/MarkDuplicates.jar
                INPUT=$input.bam
                REMOVE_DUPLICATES=true
                VALIDATION_STRINGENCY=LENIENT
                AS=true
                METRICS_FILE=$output.metrics
                CREATE_INDEX=true
                OUTPUT=$output.bam
        """
    }
}

count_reads_DNA = {
    output.dir = "counts"

    produce("counts.txt") {
        exec """
            featureCounts -O -g GeneID -T $NTHREADS -F SAF -a $SAF -o $output $inputs.bam
        """
    }
}

