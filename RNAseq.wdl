version 1.0

workflow RNAseq {
  input {
    String Queue
    String PU
    File RefFlat
    File ReferenceDir
    File? AlignedReadsReference
    File Reference
    String PL = "illumina"
    File? Fastq1
    String LB
    File OutputDir
    File? AlignedReads
    String JobGroup
    File Annotation
    String ID
    String SM
    File rRNA_intervals
    File? Fastq2
    String Name
  }
  
  if (defined(AlignedReads)) {
    call revert_aligned_reads {
      input:
        Bam = AlignedReads,
        Reference = select_first([AlignedReadsReference, Reference]),
        queue = Queue,
        jobGroup = JobGroup
    }
  }
  
  call trim_fastqs {
    input:
      F1 = select_first([revert_aligned_reads.read1, Fastq1]),
      F2 = select_first([revert_aligned_reads.read2, Fastq2]),
      queue = Queue,
      jobGroup = JobGroup
  }
  
  call gather_files {
    input:
      OutputFiles = [star_align.bam_file, star_align.bam_index, star_align.gene_counts, star_align.junction_counts, star_align.final_log, run_stringtie.stringtie_gtf, run_stringtie.stringtie_exprs, run_feature_counts.count_file, run_feature_counts.count_summary, rnaseq_metrics.metrics_out, trim_fastqs.report],
      OutputDir = OutputDir,
      queue = Queue,
      jobGroup = JobGroup
  }
  
  call rnaseq_metrics {
    input:
      in = star_align.bam_file,
      label = Name,
      RiboIntervals = rRNA_intervals,
      RefFlat = RefFlat,
      queue = Queue,
      jobGroup = JobGroup
  }
  
  call run_feature_counts {
    input:
      Bam = star_align.bam_file,
      Gtf = Annotation,
      Name = Name,
      queue = Queue,
      jobGroup = JobGroup
  }
  
  call star_align {
    input:
      Read1 = trim_fastqs.read1,
      Read2 = trim_fastqs.read2,
      ReadGroup = 'ID:' + ID + '\\tLB:' + LB + '\\tSM:' + SM + '\\tPU:' + PU + '\\tPL:' + PL,
      Name = Name,
      GenomeDir = ReferenceDir,
      Annotation = Annotation,
      queue = Queue,
      jobGroup = JobGroup
  }
  
  call run_stringtie {
    input:
      in = star_align.bam_file,
      label = Name,
      Annot = Annotation,
      queue = Queue,
      jobGroup = JobGroup
  }
  
  call remove_files as rm_files {
    input:
      files = select_all([revert_aligned_reads.read1, revert_aligned_reads.read2, trim_fastqs.read1, trim_fastqs.read2]),
      order_by = star_align.bam_file,
      queue = Queue,
      jobGroup = JobGroup
  }

}
task run_stringtie {
  input {
    String in
    String label
    String Annot
    String queue
    String jobGroup
  }

  output {
    File stringtie_gtf = "${label}.stringtie.gtf"
    File stringtie_exprs = "${label}.exprs.txt"
  }
  command <<<

    /usr/bin/stringtie ~{in} -p 4 -l ~{label} -e -G ~{Annot} -A "~{label}.exprs.txt" -o "~{label}.stringtie.gtf"
  
  >>>
  runtime {
    queue: queue
    docker_image: "mgibio/rnaseq"
    cpu: "4"
    job_group: jobGroup
    memory: "32 G"
  }

}
task gather_files {
  input {
    Array[String] OutputFiles
    String OutputDir
    String queue
    String jobGroup
  }


  output {
    String out = stdout()
  }
  command <<<

    /bin/mv -f -t ~{OutputDir} ~{sep=" "  OutputFiles}
  
  >>>
  runtime {
    queue: queue
    docker_image: "ubuntu:xenial"
    cpu: "1"
    job_group: jobGroup
    memory: "10 G"
  }

}
task run_feature_counts {
  input {
    String Bam
    String Gtf
    String Name
    String queue
    String jobGroup
  }


  output {
    File count_file = "${Name}.exon_counts.txt"
    File count_summary = "${Name}.exon_counts.txt.summary"
  }
  command <<<

    /usr/local/bin/featureCounts -T 4 -s 2 -p -B -C -t exon -g gene_id -a ~{Gtf} -o "~{Name}.exon_counts.txt" ~{Bam}
  
  >>>
  runtime {
    queue: queue
    docker_image: "dhspence/docker-subread:1.0"
    cpu: "4"
    job_group: jobGroup
    memory: "32 G"
  }

}
task revert_aligned_reads {
  input {
    File? Bam
    String Reference
    String queue
    String jobGroup
  }


  output {
    File read1 = "R1.fastq.gz"
    File read2 = "R2.fastq.gz"
  }
  command <<<

    set -eo pipefail && \
    (set -eo pipefail && /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar RevertSam INPUT=~{Bam} OUTPUT=/dev/stdout SORT_ORDER=queryname VALIDATION_STRINGENCY=SILENT | \
    /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar SamToFastq I=/dev/stdin F="R1.fastq.gz" F2="R2.fastq.gz" NON_PF=true VALIDATION_STRINGENCY=SILENT)
  
  >>>
  runtime {
    queue: queue
    docker_image: "dhspence/docker-trimgalore"
    cpu: "1"
    job_group: jobGroup
    memory: "32 G"
  }

}
task star_align {
  input {
    File Read1
    File Read2
    String ReadGroup
    String Name
    String GenomeDir
    String Annotation
    String queue
    String jobGroup
  }


  output {
    File bam_file = "${Name}.Aligned.sortedByCoord.out.bam"
    File bam_index = "${Name}.Aligned.sortedByCoord.out.bam.bai"
    File gene_counts = "${Name}.ReadsPerGene.out.tab"
    File junction_counts = "${Name}.SJ.out.tab"
    File final_log = "${Name}.Log.final.out"
  }
  command <<<

    /usr/local/STAR/bin/Linux_x86_64/STAR  --genomeDir ~{GenomeDir} --readFilesIn ~{sep=","  Read1} ~{sep=","  Read2} \
    --outSAMunmapped Within --outSAMmapqUnique 60 --outSAMattributes NH HI AS NM MD --outSAMattrIHstart 0 --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFileNamePrefix ./"~{Name}." \
    --chimSegmentMin 12 --chimSegmentReadGapMax parameter 3 --alignSJstitchMismatchNmax 5 -1 5 5 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    --sjdbGTFfile ~{Annotation} --outTmpDir /tmp/~{Name} --sjdbFileChrStartEnd - --limitSjdbInsertNsj 3000000 --sjdbInsertSave Basic \
    --outSAMattrRGline "~{ReadGroup}" \
    --readFilesCommand zcat  --outFilterMismatchNoverReadLmax 0.05 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --sjdbOverhang 149 --outFilterScoreMinOverLread 0.33 \
    --outFilterMatchNminOverLread 0.33 --alignSoftClipAtReferenceEnds No --chimJunctionOverhangMin 15  --runThreadN 8 --genomeLoad NoSharedMemory --limitBAMsortRAM 100000000000 \
    --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outSAMheaderHD @HD VN:1.4 SO:coordinate --quantMode GeneCounts && \
    /usr/local/bin/samtools index *.bam
  
  >>>
  runtime {
    queue: queue
    docker_image: "dhspence/docker-star"
    cpu: "8"
    job_group: jobGroup
    memory: "84 G"
  }

}
task trim_fastqs {
  input {
    File F1
    File F2
    String queue
    String jobGroup
  }


  output {
    File report = glob("*.html")[0]
    File read1 = glob("*R1*trim.fastq.gz")[0]
    File read2 = glob("*R2*trim.fastq.gz")[0]
  }
  command <<<

    fastp -w 4 -i ~{F1} -o $(basename ~{F1})".trim.fastq.gz" -I ~{F2} -O $(basename ~{F2})".trim.fastq.gz"
  
  >>>
  runtime {
    queue: queue
    docker_image: "dhspence/docker-fastp"
    cpu: "4"
    job_group: jobGroup
    memory: "16 G"
  }

}
task rnaseq_metrics {
  input {
    String in
    String label
    String RiboIntervals
    String RefFlat
    String queue
    String jobGroup
  }


  output {
    File metrics_out = "${label}.rnaseq_metrics.txt"
  }
  command <<<

    /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectRnaSeqMetrics I=~{in} O="~{label}.rnaseq_metrics.txt" STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
    REF_FLAT=~{RefFlat} RIBOSOMAL_INTERVALS=~{RiboIntervals} VALIDATION_STRINGENCY=LENIENT
  
  >>>
  runtime {
    queue: queue
    docker_image: "registry.gsc.wustl.edu/genome/picard-2.12.1-r:1"
    cpu: "1"
    job_group: jobGroup
    memory: "10 G"
  }

}
task remove_files {
  input {
    Array[String] files
    String order_by
    String queue
    String jobGroup
  }


  output {
    String out = stdout()
  }
  command <<<

    /bin/rm -f ~{sep=" "  files}
  
  >>>
  runtime {
    queue: queue
    docker_image: "ubuntu:xenial"
    cpu: "1"
    job_group: jobGroup
    memory: "10 G"
  }

}

