fastqc = {
    doc "Quality control using FASTQC"
    output.dir = "fastqc"

    transform(".fastq.gz") to ("_fastqc.zip") {
        exec "fastqc -o $output.dir $input.gz"
    }
}

trim = {
    doc "Trim poor quility bases and/or adapter sequences from reads using Trimmomatic"
    output.dir = "trimmed"

    filter("trim") {
        exec """ 
            trimmomatic SE -phred33 $input.gz $output.gz ILLUMINACLIP:$ADAPTERS:2:30:10 
            TRAILING:20 LEADING:20 MINLEN:35
        """
    }
}

trim_PE = {
   doc "Trim reads using Trimmomatic"
   output.dir="trimmed"
   filter("trim","trim") {
        exec """
                trimmomatic PE -phred33
                $input1.gz $input2.gz
                $output1.gz ${output1.prefix}.unpaired.gz
                $output2.gz ${output2.prefix}.unpaired.gz
                LEADING:15 TRAILING:15 SLIDINGWINDOW:5:15 MINLEN:30
        """
       }
}

sort_bam = {
    doc "Sort bam files"
    output.dir = "sorted_indexed"

    filter("srt") {
        exec "samtools sort $input.bam $output.prefix"
    }
}

index_bam = {
    doc "Index bam files"
    output.dir = "sorted_indexed"

    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam"
    }
    forward input
}

