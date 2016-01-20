star_map_1pass_SE = {
    String files = "${inputs}";
    files = files.replaceAll(' ',',');

    doc "Map reads using the STAR aligner: 1st pass"
    output.dir = MAPPING_DIR+"/1pass"

    from("*.fastq.trim.gz") {
        produce ("SJ.out.tab") {
            exec """
                STAR --genomeDir $GEN_DIR --readFilesIn ${files} --readFilesCommand zcat --outSAMtype None
                --runThreadN $NTHREADS --outFileNamePrefix ${output.dir}/
            """
        }
    }
}

star_gen2pass = {
    doc "Map reads using the STAR aligner: generate genome"
    output.dir = MAPPING_DIR+"/genome_2pass"

    produce("Genome"){
        exec """
            STAR --runMode genomeGenerate --genomeDir ${output.dir} --genomeFastaFiles $GEN_FASTA --runThreadN $NTHREADS
            --sjdbFileChrStartEnd ${MAPPING_DIR}/1pass/SJ.out.tab --sjdbOverhang 99 --sjdbGTFfile $GTF
            --outFileNamePrefix ${output.dir}/
        """
    }
}

star_map_2pass_SE = {
    String files = "${inputs}";
    files = files.replaceAll(' ',',');

    doc "Map reads using the STAR aligner: 2nd pass"
    output.dir = MAPPING_DIR

    // ensure that filename suffixes in transform statement match *actual* filenames
    transform("(.*)_L001_R1_001.fastq.trim.gz","(.*)_L002_R1_001.fastq.trim.gz",
		"(.*)_L003_R1_001.fastq.trim.gz","(.*)_L004_R1_001.fastq.trim.gz") to ("\$1.Aligned.sortedByCoord.out.bam") {
        exec """
	      STAR --genomeDir mapped/genome_2pass --readFilesIn ${files} \
	      --outFileNamePrefix ${output.prefix.prefix.prefix.prefix}.   \
	      --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --runThreadN 5
      	"""
   }
}


count_reads_RNA = {
    doc "Count reads"
    output.dir = "counts"

    produce("counts.txt") {
        exec """
            featureCounts --primary -p -t exon -g GeneID -T $NTHREADS -F SAF -a $SAF -o $output $inputs.bam
        """
    }
}


