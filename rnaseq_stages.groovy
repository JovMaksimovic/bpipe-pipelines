star_genome_gen = {
    //Generate STAR genome index
    doc "Generate STAR genome index"

    produce("Genome") {
            exec """
                $STAR --runMode genomeGenerate 
		--runThreadN $NTHREADS
		--genomeDir $GEN_DIR
		--genomeFastaFiles $GEN_FASTA
		--sjdbGTFfile $GTF
		--sjdbOverhang $SJBOHANG
            ""","stargenind"
    }
}

star_map_1pass_SE = {
    //Map single-end reads using the STAR aligner: 1st pass
    String files = "${inputs}";
    files = files.replaceAll(' ',',');

    doc "Map single-end reads using the STAR aligner: 1st pass"
    output.dir = "mapped/1pass"

    from("*.fastq.gz") {
        produce ("SJ.out.tab") {
            exec """
                $STAR --genomeDir $GEN_DIR --readFilesIn ${files} --readFilesCommand zcat --outSAMtype None
                --runThreadN $NTHREADS --sjdbGTFfile $GTF --outFileNamePrefix ${output.dir}/
            ""","star1pass"
        }
    }
}

star_map_1pass_PE = {
  //Map paired-end reads using the STAR aligner: 1st pass
  println inputs;
  String R1 = inputs[0];
  String R2 = inputs[1];

  for (int i = 2; i < inputs.size(); i++) {
        if (i % 2 == 0){
                R1 = R1+","+inputs[i];
        } else {
                R2 = R2+","+inputs[i];
        }
  }
  
  String files = R1+" "+R2;
  println files;

  doc "Map paired-end reads using the STAR aligner: 1st pass"
  output.dir = "mapped/1pass"

  from("*.fastq.trim.gz") {
    produce ("SJ.out.tab") {
      exec """
        $STAR --genomeDir $GEN_DIR --readFilesIn ${files} --readFilesCommand zcat --outSAMtype None
        --runThreadN $NTHREADS --sjdbGTFfile $GTF --outFileNamePrefix ${output.dir}/
      ""","star1pass"
    }
  }
}

//star_gen2pass = {
//    doc "Map reads using the STAR aligner: generate genome"
//    output.dir = "mapped/genome_2pass"

//    produce("Genome"){
//        exec """
//            	$STAR --runMode genomeGenerate --genomeDir ${output.dir} --genomeFastaFiles $GEN_FASTA 
//		--runThreadN $NTHREADS --sjdbFileChrStartEnd ${output.dir}/1pass/SJ.out.tab 
//		--sjdbOverhang $SJBOHANG --sjdbGTFfile $GTF --outFileNamePrefix ${output.dir}/
 //       ""","stargen"
//    }
//}

star_map_2pass_SE = {
    //Map single-end reads using the STAR aligner; 2nd pass
    String files = "${inputs}";
    files = files.replaceAll(' ',',');
    
    doc "Map single-end reads using the STAR aligner: 2nd pass"
    output.dir = "mapped"

    transform("trim.fastq.gz") to ("Aligned.sortedByCoord.out.bam") {
        exec """
	      	$STAR --genomeDir $GEN_DIR --readFilesIn ${files} 
	      	--outFileNamePrefix ${output.prefix.prefix.prefix}.   
		--sjdbFileChrStartEnd ${output.dir}/1pass/SJ.out.tab --sjdbOverhang $SJBOHANG
	      	--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --runThreadN $NTHREADS
      	""","star2pass"
   }
}

star_map_2pass_PE = {
  //Map paired-end reads using the STAR aligner: 2nd pass
  doc "Map paired-end reads using the STAR aligner: 2nd pass"
  output.dir = "mapped"

  transform("(.*)_1.fastq.trim.gz","(.*)_2.fastq.trim.gz") to ("\$1.Aligned.out.bam") {
    exec """
        $STAR --genomeDir $GEN_DIR --readFilesIn $input1.gz $input2.gz
	--sjdbFileChrStartEnd ${output.dir}/1pass/SJ.out.tab --sjdbOverhang $SJBOHANG
        --outFileNamePrefix ${output.prefix.prefix.prefix}. --readFilesCommand zcat
        --outSAMtype BAM Unsorted --runThreadN $NTHREADS
    ""","star2pass"
  }
}

count_reads_RNA = {
    //Count reads across features using RNA data
    doc "Count reads across features from RNA data"
    output.dir = "counts"

    produce("counts.txt") {
        exec """
            $FEATURECOUNTS --primary -p -t exon -g GeneID -T $NTHREADS -F SAF -a $SAF -o $output $inputs.bam
        ""","count"
    }
}


